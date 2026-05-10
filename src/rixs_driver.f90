!> RIXS driver: compute the Resonant Inelastic X-ray Scattering spectrum.
!!
!! Called by rixs_fsolver (the f2py-exposed wrapper in pyapi.f90), which is in
!! turn invoked from the Python solvers rixs_siam_fort, rixs_1v1c_fort and
!! rixs_2v1c_fort defined in edrixs/solvers.py.  The Python layer loops over
!! incident photon energies omega_in and over polarisation channels, rewriting
!! config.in (with omega_in, gamma_in) and the transition operator files
!! transop_rixs_i.in / transop_rixs_f.in for each combination.  Hence
!! rixs_driver is called once per (omega_in, polarisation) pair.
!!
!! Implements the Kramers-Heisenberg RIXS cross-section of equation (2) of
!! Wang et al. (CPC 2019),
!!   I_RIXS(omega_in, omega, k_i, k_f, eps_i, eps_f)
!!     = sum_i (1/Z) e^{-E_i / K_B T} sum_f
!!         |<f| hat{D}_f^dagger (omega_in - H_n + E_i + i*Gamma_c)^{-1} hat{D}_i |i>|^2
!!         * Gamma/pi / [(omega - E_f + E_i)^2 + Gamma^2],
!! using the algorithm of section 2.4 of the paper.  The sum over final
!! states f is reformulated as the imaginary part of a matrix element of a
!! resolvent in the final-state space, then evaluated as a continued fraction
!! via Lanczos (equation (12) of the paper):
!!   I_RIXS(omega_in, omega) = <F| (omega - H_i + E_i + i*Gamma)^{-1} |F>.
!!
!! Workflow (executed once per ground state igs = 1..num_gs, mirroring the
!! flow diagram in Fig. 2 of the paper):
!!  1. Build the photon-absorption transition operator hat{D}_i as a distributed
!!     CSR matrix (ndim_n x ndim_i) using build_transop_i, then apply it to the
!!     i-th ground state |Gamma_i> (read from eigvec.igs):
!!         |b> = hat{D}_i |Gamma_i>      (paper equation following (10))
!!  2. Build the intermediate-state Hamiltonian H_n with build_ham_n, but with
!!     ham_sign = -1.0 and a complex shift omega = omega_in + E_i + i*Gamma_c.
!!     This stores the matrix -(H_n - omega I), so pminres_csr (equation (11)
!!     of the paper) solves
!!         (omega_in - H_n + E_i + i*Gamma_c) |x> = |b>.
!!  3. Solve the linear system with pminres_csr to obtain
!!         |x> = (omega_in - H_n + E_i + i*Gamma_c)^{-1} |b>.
!!     This is the most expensive step of the RIXS calculation; convergence
!!     is controlled by linsys_tol and linsys_max.
!!  4. Build the photon-emission transition operator hat{D}_f^dagger (ndim_f x ndim_n)
!!     and apply it to |x>:  |F> = hat{D}_f^dagger |x>  (paper equation following (11)).
!!     This is the seed vector for the final-state Krylov expansion.
!!  5. Build the final-state Hamiltonian H_f.  Because final states are
!!     eigenstates of H_i (no core hole), H_f reuses build_ham_i with the
!!     hopping_i / coulomb_i parameters.
!!  6. Run build_krylov_mp on H_i starting from |F> for up to nkryl Lanczos
!!     steps.  Save alpha, beta and ||F||^2 to rixs_poles.igs; Python then
!!     evaluates the paper's continued-fraction G(omega) = <F|(omega - H_i +
!!     E_i + i*Gamma)^{-1}|F> on the energy-loss mesh.
!!
!! Input files (written by Python before rixs_fsolver is called):
!!  - config.in            : control namelist with omega_in, gamma_in, ...
!!  - hopping_i.in, coulomb_i.in : initial/final-state Hamiltonian (no core hole)
!!  - hopping_n.in, coulomb_n.in : intermediate-state Hamiltonian (with core hole)
!!  - transop_rixs_i.in    : absorption transition operator (initial -> intermediate)
!!  - transop_rixs_f.in    : emission transition operator (intermediate -> final)
!!  - fock_i.in, fock_n.in, fock_f.in : Fock bases for the three sectors
!!  - eigvec.igs           : initial-state eigenvector(s) from a prior ED run
!!
!! Output files (read by Python):
!!  - rixs_poles.igs       : Lanczos coefficients alpha, beta, norm and E_gs;
!!                            evaluated by build_spectrum/get_spectra_from_poles
!!                            in Python to give I(omega_in, omega_loss).
subroutine rixs_driver()
    use m_constants
    use m_control
    use m_types
    use m_global
    use m_lanczos
    use mpi

    implicit none

    integer :: nblock
    integer :: mloc
    integer :: nloc
    integer :: info
    integer :: neff
    integer :: igs
    integer :: needed(nprocs,nprocs)
    integer :: needed2(nprocs,nprocs)
    integer :: end_indx(2,2,nprocs)
    integer :: end_indx2(2,2,nprocs)
    integer :: ierror
    integer(dp) :: num_of_nonzeros

    real(dp) :: rtemp
    real(dp) :: norm
    real(dp) :: eigvals

    complex(dp)              :: omega
    complex(dp), allocatable :: eigvecs(:)
    complex(dp), allocatable :: eigvecs_mpi(:)
    complex(dp), allocatable :: phi_vec(:)
    complex(dp), allocatable :: x_vec(:)

    character(len=20) :: fname
    character(len=10) :: char_I

    call read_hopping_i()
    call read_coulomb_i()
    call read_hopping_n()
    call read_coulomb_n()
    call read_transop_rixs_i()
    call read_transop_rixs_f()
    call read_fock_i()
    call read_fock_n()
    call read_fock_f()

    ndim_n = ndim_n_nocore * num_core_orbs

    if (myid == master) then
        print *, " fedrixs >>> RIXS Begin ..."
        print *
        write(mystd,"(a20,i15)")    "num_val_orbs:  ", num_val_orbs
        write(mystd,"(a20,i15)")    "num_core_orbs: ", num_core_orbs
        write(mystd,"(a20,i15)")    "ndim_i:        ", ndim_i
        write(mystd,"(a20,i15)")    "ndim_n:        ", ndim_n
        write(mystd,"(a20,i15)")    "ndim_f:        ", ndim_f
        write(mystd,"(a20,i15)")    "nhopp_i:       ", nhopp_i
        write(mystd,"(a20,i15)")    "nhopp_n:       ", nhopp_n
        write(mystd,"(a20,i15)")    "ncoul_i:       ", ncoul_i
        write(mystd,"(a20,i15)")    "ncoul_n:       ", ncoul_n
        write(mystd,"(a20,i15)")    "num_gs:        ", num_gs
        write(mystd,"(a20,i15)")    "nkryl:         ", nkryl
        write(mystd,"(a20,i15)")    "linsys_max:    ", linsys_max
        write(mystd,"(a20,e15.2)")  "linsys_tol:    ", linsys_tol
        write(mystd,"(a20,f15.6)")  "omega_in:      ", omega_in
        write(mystd,"(a20,f15.6)")  "gamma_in       ", gamma_in
        print *
    endif

    call dealloc_fock_i()
    call dealloc_fock_n()
    call dealloc_fock_f()
    do igs=1, num_gs
        if (myid == master) then
            print *, " fedrixs >>> For initial state:  ", igs
        endif
        if (myid == master) then
            print *, "    Building transition operator for absorption process ..."
        endif

        nblock = nprocs
        call partition_task(nprocs, ndim_n, ndim_i, end_indx)
        mloc = end_indx(2,1,myid+1)-end_indx(1,1,myid+1) + 1
        nloc = end_indx(2,2,myid+1)-end_indx(1,2,myid+1) + 1
        call read_fock_i()
        call read_fock_n()
        call alloc_tran_csr(nblock)
        call build_transop_i(ndim_n_nocore, ndim_i, fock_n, fock_i, num_val_orbs,&
          num_core_orbs, nblock, end_indx, ntran_rixs_i, transop_rixs_i, tran_csr)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_fock_i()
        call dealloc_fock_n()
        call get_needed_indx(nblock, tran_csr, needed)
        if (myid==master) then
            print *, "    Done !"
            print *
        endif

        allocate(eigvecs(nloc))
        allocate(eigvecs_mpi(ndim_i))
        eigvecs_mpi = czero
        eigvecs = czero
        ! read eigenvectors from files generated by ed.x
        write(char_I, '(i5)') igs
        fname="eigvec."//trim(adjustl(char_I))
        call read_eigvecs(fname, ndim_i, eigvecs_mpi, eigvals)
        eigvecs = eigvecs_mpi(end_indx(1,2,myid+1): end_indx(2,2,myid+1))
        deallocate(eigvecs_mpi)
        if (myid==master) then
            print *, "    Apply transition operator on the ground state to get intermediate state..."
        endif

        allocate(phi_vec(mloc))
        phi_vec = czero
        call pspmv_csr(new_comm, nblock, end_indx, needed, mloc, nloc, tran_csr, eigvecs, phi_vec)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_tran_csr(nblock)
        deallocate(eigvecs)
        if (myid==master) then
            print *, "    Done !"
            print *
        endif
        ! build the Hamiltonian for intermediate states
        if (myid==master) then
            print *, "    Building Hamiltonian for intermediate configuration ..."
        endif
        call partition_task(nprocs, ndim_n, ndim_n, end_indx2)
        call read_fock_n()
        call alloc_ham_csr(nblock)
        rtemp = -1.0_dp
        omega = dcmplx(omega_in+eigvals, gamma_in)
        call build_ham_n(ndim_n_nocore, fock_n, num_val_orbs, num_core_orbs, nblock, end_indx2, &
                        nhopp_n, hopping_n, ncoul_n, coulomb_n, omega, rtemp, ham_csr)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_fock_n()
        call get_needed_indx(nblock, ham_csr, needed2)
        call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
        if (myid==master) then
            print *, "    Number of nonzero elements of intermediate Hamiltonian: ", num_of_nonzeros
            print *, "    Done !"
            print *
        endif

        ! solve linear equation
        if (myid==master) then
            print *, "    Solve linear equation ..."
        endif
        allocate(x_vec(mloc))
        x_vec = czero
        call pminres_csr(nblock, end_indx2, needed2, mloc, ham_csr, phi_vec, x_vec, info)
        call MPI_BARRIER(new_comm, ierror)
        deallocate(phi_vec)
        call dealloc_ham_csr(nblock)
        if (myid==master) then
            print *, "    Done !"
            print *
        endif

        ! get the final state
        if (myid==master) then
            print *, "    Building transition operator for the emission process ..."
        endif
        call partition_task(nprocs, ndim_f, ndim_n, end_indx)
        mloc = end_indx(2,1,myid+1)-end_indx(1,1,myid+1) + 1
        nloc = end_indx(2,2,myid+1)-end_indx(1,2,myid+1) + 1
        call alloc_tran_csr(nblock)
        call read_fock_f()
        call read_fock_n()
        call build_transop_f(ndim_f, ndim_n_nocore, fock_f, fock_n, num_val_orbs, &
           num_core_orbs, nblock, end_indx, ntran_rixs_f, transop_rixs_f, tran_csr)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_fock_f()
        call dealloc_fock_n()
        call get_needed_indx(nblock, tran_csr, needed)
        if (myid==master) then
            print *, "    Done !"
            print *
        endif

        if (myid==master) then
            print *, "    Applying transition operator to get final state ..."
        endif
        allocate(phi_vec(mloc))
        phi_vec = czero
        call pspmv_csr(new_comm, nblock, end_indx, needed, mloc, nloc, tran_csr, x_vec, phi_vec)
        call MPI_BARRIER(new_comm, ierror)
        deallocate(x_vec)
        call dealloc_tran_csr(nblock)

        ! build the Hamiltonian for final states
        if (myid==master) then
            print *, "    Building Hamiltonian for initial configuration ..."
        endif
        call partition_task(nprocs, ndim_f, ndim_f, end_indx2)
        call alloc_ham_csr(nblock)
        rtemp = 1.0_dp
        omega = dcmplx(0.0_dp, 0.0_dp)
        call read_fock_f()
        call build_ham_i(ndim_f, fock_f, nblock, end_indx2, nhopp_i, hopping_i, ncoul_i, coulomb_i, omega, rtemp, ham_csr)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_fock_f()
        call get_needed_indx(nblock, ham_csr, needed2)
        call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
        if (myid==master) then
            print *, "    Number of nonzero elements of initial Hamiltonian: ", num_of_nonzeros
            print *, "    Done !"
            print *
        endif

        ! get the rixs spectrum
        if (myid==master) then
            print *, "    Building Krylov space for RIXS spectrum ..."
        endif
        allocate(krylov_alpha(nkryl))
        allocate(krylov_beta(0:nkryl))
        krylov_alpha = zero
        krylov_beta  = zero
        call build_krylov_mp(nblock, end_indx2, needed2, mloc, ham_csr, phi_vec, nkryl, neff, krylov_alpha, krylov_beta, norm)
        call MPI_BARRIER(new_comm, ierror)
        call dealloc_ham_csr(nblock)
        deallocate(phi_vec)
        write(char_I, '(I5)') igs
        fname="rixs_poles."//trim(adjustl(char_I))
        call write_krylov(fname, neff, krylov_alpha(1:neff), krylov_beta(1:neff), norm, eigvals)
        deallocate(krylov_alpha)
        deallocate(krylov_beta)
        if (myid==master) then
            print *, "    Done !"
            print *
        endif
    enddo ! igs=1, num_gs

    call dealloc_transop_rixs_i()
    call dealloc_transop_rixs_f()
    call dealloc_hopping_i()
    call dealloc_hopping_n()
    call dealloc_coulomb_i()
    call dealloc_coulomb_n()

    if (myid==master) then
        print *
        print *, " fedrixs >>> RIXS End !"
    endif

    return
end subroutine rixs_driver
