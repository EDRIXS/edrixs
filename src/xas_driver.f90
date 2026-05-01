!> XAS driver: compute the X-ray Absorption Spectrum via the Krylov method.
!!
!! For each of the num_gs initial states (stored as eigvec.k files from ed.x):
!!
!!  1. Build the absorption transition operator T as a distributed CSR matrix
!!     mapping the initial-state Hilbert space into the intermediate-state space
!!     (ndim_n rows, ndim_i columns).
!!  2. Apply T to the ground-state eigenvector to produce phi = T|psi_0>, the
!!     seed vector in the intermediate-state space.
!!  3. Build and diagonalise the intermediate-state Hamiltonian H_n.
!!  4. Run build_krylov_mp to compute the Lanczos coefficients alpha/beta that
!!     define the continued-fraction Green's function G(omega) = <phi|(omega-H_n)^{-1}|phi>.
!!  5. Write the Krylov data to xas_poles.k for post-processing with build_spectrum.
!!
!! The XAS intensity for a given polarisation is proportional to -Im G(omega_in).
subroutine xas_driver()
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
    complex(dp), allocatable :: phi_vec(:)
    complex(dp), allocatable :: eigvecs_mpi(:)

    character(len=20) :: fname
    character(len=10) :: char_I

    if (myid == master) then
        print *, "---------------------------"
        print *, " fedrixs >>> XAS Begin ... "
        print *
    endif

    call read_hopping_n()
    call read_coulomb_n()
    call read_transop_xas()
    call read_fock_i()
    call read_fock_n()

    ndim_n = ndim_n_nocore * num_core_orbs

    if (myid == master) then
       print *
       write(mystd,"(a20, i15)")  " num_val_orbs:  ", num_val_orbs
       write(mystd,"(a20, i15)")  " num_core_orbs: ", num_core_orbs
       write(mystd,"(a20, i15)")  " ndim_i:        ", ndim_i
       write(mystd,"(a20, i15)")  " ndim_n:        ", ndim_n
       write(mystd,"(a20, i15)")  " nhopp_n:       ", nhopp_n
       write(mystd,"(a20, i15)")  " ncoul_n:       ", ncoul_n
       write(mystd,"(a20, i15)")  " num_gs:        ", num_gs
       write(mystd,"(a20, i15)")  " nkryl:         ", nkryl
       print *
    endif
    call dealloc_fock_i()
    call dealloc_fock_n()

    do igs=1, num_gs
        if (myid == master) then
            print *, " fedrixs >>> For ground state: ", igs
            print *, "    Building transition operator for absorption process ..."
        endif
        call read_fock_i()
        call read_fock_n()
        nblock = nprocs
        call partition_task(nprocs, ndim_n, ndim_i, end_indx)
        mloc = end_indx(2,1,myid+1)-end_indx(1,1,myid+1) + 1
        nloc = end_indx(2,2,myid+1)-end_indx(1,2,myid+1) + 1
        call alloc_tran_csr(nblock)
        call build_transop_i(ndim_n_nocore, ndim_i, fock_n, fock_i, num_val_orbs, &
                 num_core_orbs, nblock, end_indx, ntran_xas, transop_xas, tran_csr)
        call MPI_BARRIER(new_comm, ierror)

        call get_needed_indx(nblock, tran_csr, needed)
        call dealloc_fock_i()
        if (myid == master) then
            print *, "    Done !"
            print *
        endif

        if (myid==master) then
            print *, "    Apply transition operator on the ground state to get intermediate state..."
        endif
        allocate(eigvecs(nloc))
        allocate(eigvecs_mpi(ndim_i))
        eigvecs = czero
        eigvecs_mpi = czero
        eigvals = zero
        write(char_I, '(i5)') igs
        fname="eigvec."//trim(adjustl(char_I))
        call read_eigvecs(fname, ndim_i, eigvecs_mpi, eigvals)
        eigvecs = eigvecs_mpi(end_indx(1,2,myid+1): end_indx(2,2,myid+1))
        deallocate(eigvecs_mpi)
        allocate(phi_vec(mloc))
        phi_vec = czero

        call MPI_BARRIER(new_comm, ierror)
        call pspmv_csr(new_comm, nblock, end_indx, needed, mloc, nloc, tran_csr, eigvecs, phi_vec)
        call MPI_BARRIER(new_comm, ierror)
        deallocate(eigvecs)
        call dealloc_tran_csr(nblock)

        if (myid==master) then
            print *, "    Done !"
            print *
        endif

        if (myid == master) then
            print *, "    Build Hamiltonian for intermediate configuration..."
        endif
        call partition_task(nprocs, ndim_n, ndim_n, end_indx2)
        call alloc_ham_csr(nblock)
        rtemp = 1.0_dp
        omega = dcmplx(0.0_dp,0.0_dp)
        call build_ham_n(ndim_n_nocore, fock_n, num_val_orbs, num_core_orbs, nblock, end_indx2, &
                         nhopp_n, hopping_n, ncoul_n, coulomb_n, omega, rtemp, ham_csr)
        call MPI_BARRIER(new_comm, ierror)
        call get_needed_indx(nblock, ham_csr, needed2)
        call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
        call dealloc_fock_n()
        if (myid == master) then
            print *, "    Number of nonzero elements of intermediate Hamiltonian: ", num_of_nonzeros
            print *, "    Done !"
            print *
        endif

        if (myid == master) then
            print *, "    Building Krylov space for XAS spectrum ..."
        endif
        allocate(krylov_alpha(nkryl))
        allocate(krylov_beta(0:nkryl))
        krylov_alpha = zero
        krylov_beta  = zero
        call build_krylov_mp(nblock, end_indx2, needed2, mloc, ham_csr, phi_vec, nkryl, neff, krylov_alpha, krylov_beta, norm)

        write(char_I, '(I5)') igs
        fname="xas_poles."//trim(adjustl(char_I))
        call write_krylov(fname, neff, krylov_alpha(1:neff), krylov_beta(1:neff), norm, eigvals)

        call dealloc_ham_csr(nblock)
        deallocate(krylov_alpha)
        deallocate(krylov_beta)
        deallocate(phi_vec)
        if (myid == master) then
            print *, "    Done !"
            print *
        endif
    enddo ! over igs=1, num_gs

    call dealloc_transop_xas()
    call dealloc_hopping_n()
    call dealloc_coulomb_n()
    if (myid == master) then
        print *
        print *, " fedrixs >>> XAS End !"
    endif

    return
end subroutine xas_driver
