!> Exact Diagonalization (ED) driver: diagonalise the initial-state Hamiltonian.
!!
!! Called by ed_fsolver (the f2py-exposed wrapper in pyapi.f90), which is in
!! turn invoked from the Python solvers ed_siam_fort, ed_1v1c_fort and
!! ed_2v1c_fort defined in edrixs/solvers.py.  The same kernel serves all three
!! Python entry points; only the Hamiltonian elements written to the input
!! files differ (atomic many-body for ed_1v1c_fort/ed_2v1c_fort, single-impurity
!! Anderson model for ed_siam_fort).
!!
!! Workflow:
!!  1. read_hopping_i, read_coulomb_i, read_fock_i load the initial Hamiltonian
!!     parameters from hopping_i.in, coulomb_i.in, fock_i.in (written by Python).
!!  2. partition_task splits the ndim_i x ndim_i Hilbert space across ranks.
!!  3. build_ham_i constructs the distributed CSR Hamiltonian.
!!  4. The eigensolver selected by ed_solver runs:
!!       - 0: full LAPACK ZHEEV (gathered to master, broadcast back)
!!       - 1: parallel Lanczos (no re-orthogonalisation, less accurate)
!!       - 2: parallel ARPACK   (recommended; default for SIAM)
!!     If ndim_i < min_ndim, ed_solver is forced to 0.
!!  5. The single-particle density matrix is computed for each of the nvector
!!     lowest eigenstates by accumulating contributions <Gamma_i|f^dagger_alpha f_beta|Gamma_i>
!!     over the local Fock basis (notation from Table 1 of Wang et al. CPC 2019).
!!     When ed_solver != 0 the eigenvectors must
!!     first be gathered with MPI_ISEND/MPI_RECV because off-diagonal density
!!     matrix elements connect Fock states owned by different ranks.
!!
!! Input files (written by Python in solvers.py):
!!  - config.in        : Fortran namelist with ed_solver, num_val_orbs, neval, ...
!!  - hopping_i.in     : non-zero one-body matrix elements
!!  - coulomb_i.in     : non-zero two-body (Coulomb) matrix elements
!!  - fock_i.in        : Fock basis (bit-string integers)
!!
!! Output files (read back by Python):
!!  - eigvals.dat      : neval lowest eigenvalues
!!  - denmat.dat       : single-particle density matrices for the nvector states
!!  - eigvec.k         : (if idump=.true.) full eigenvector for state k,
!!                       consumed later by xas_driver/rixs_driver
subroutine ed_driver()
    use m_constants
    use m_control
    use m_types
    use m_global
    use m_lanczos
    use mpi

    implicit none

    integer :: i,j,k
    integer :: icfg, jcfg
    integer :: mblock
    integer :: nblock
    integer :: nloc
    integer :: needed(nprocs,nprocs)
    integer :: end_indx(2,2,nprocs)
    integer :: ierror
    integer :: info
    integer :: actual_step
    integer :: nconv
    integer :: request
    integer :: stat(MPI_STATUS_SIZE)
    integer :: sgn
    integer :: tot_sign
    integer(dp) :: old, new
    integer(dp) :: num_of_nonzeros
    integer, external :: binary_search

    real(dp) :: rtemp
    real(dp) :: time_begin
    real(dp) :: time_end
    real(dp) :: time_used
    real(dp), allocatable :: eigvals(:)
    real(dp), allocatable :: eigval_full(:)

    complex(dp) :: omega
    complex(dp) :: denmat(num_val_orbs, num_val_orbs, nvector)
    complex(dp) :: denmat_mpi(num_val_orbs, num_val_orbs, nvector)

    complex(dp), allocatable :: ham_full(:,:)
    complex(dp), allocatable :: ham_full_mpi(:,:)
    complex(dp), allocatable :: eigvecs_full(:,:)

    complex(dp), allocatable :: eigvecs(:,:)
    complex(dp), allocatable :: eigvecs_mpi(:)

    character(len=20) :: fname
    character(len=10) :: char_I

    time_begin = 0.0_dp
    time_end   = 0.0_dp
    time_used  = 0.0_dp
    call cpu_time(time_begin)

    if (myid == master) then
        print *, "--------------------------------------------"
        print *, " fedrixs >>> ED begin ... "
        print *
    endif
    call read_hopping_i()
    call read_coulomb_i()
    call read_fock_i()

    if (ndim_i < min_ndim ) then
        if (myid==master) then
            print *, " fedrixs >>> ndim_i:", ndim_i, " is smaller than min_ndim:", min_ndim
            print *, " fedrixs >>> set ed_solver = 0, use full-diagonalization !"
            print *
        endif
        ed_solver = 0
    endif

    if(ndim_i < neval) then
        neval = ndim_i
    endif
    if(ndim_i < nvector) then
        nvector=ndim_i
    endif

    if (myid == master) then
        write(mystd,"(a20, i15)")   "num_val_orbs:  ", num_val_orbs
        write(mystd,"(a20, i15)")   "ed_solver:     ", ed_solver
        write(mystd,"(a20, i15)")   "neval:         ", neval
        write(mystd,"(a20, i15)")   "nvector:       ", nvector
        write(mystd,"(a20, i15)")   "maxiter:       ", maxiter
        write(mystd,"(a20, i15)")   "min_ndim:      ", min_ndim
        write(mystd,"(a20, i15)")   "ncv:           ", ncv
        write(mystd,"(a20, i15)")   "nhopp_i:       ", nhopp_i
        write(mystd,"(a20, i15)")   "ncoul_i:       ", ncoul_i
        write(mystd,"(a20, i15)")   "ndim_i:        ", ndim_i
        write(mystd,"(a20, i15)")   "nprocs:        ", nprocs

        write(mystd,"(a20, e15.2)") "eigval_tol:    ", eigval_tol
        write(mystd,"(a20, l15)")   "idump:         ", idump
        print *
    endif

    if (myid==master) then
        print *, " fedrixs >>> Build Hamiltonian ..."
    endif
    call partition_task(nprocs, ndim_i, ndim_i, end_indx)
    rtemp = 1.0_dp
    omega = czero
    mblock = nprocs
    nblock = nprocs
    nloc = end_indx(2,2,myid+1)-end_indx(1,2,myid+1) + 1
    call alloc_ham_csr(nblock)
    call build_ham_i(ndim_i, fock_i, nblock, end_indx, nhopp_i, hopping_i, ncoul_i, coulomb_i, omega, rtemp, ham_csr)
    call MPI_BARRIER(new_comm, ierror)
    call get_needed_indx(nblock, ham_csr, needed)
    call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
    call dealloc_fock_i()
    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Number of nonzero elements of the Hamiltonian", num_of_nonzeros
        print *, " fedrixs >>> Done ! Time used: ", time_end - time_begin, "  seconds"
        print *
        print *, " fedrixs >>> Diagonalize Hamiltonian to find a few lowest states ..."
        print *
    endif
    time_begin = time_end

    allocate(eigvals(neval))
    ! full diagonalization
    if (ed_solver == 0) then
        allocate(ham_full(ndim_i, ndim_i))
        allocate(ham_full_mpi(ndim_i, ndim_i))
        allocate(eigvecs_full(ndim_i, ndim_i))
        allocate(eigval_full(ndim_i))
        allocate(eigvecs(ndim_i,nvector))
        ham_full = czero
        ham_full_mpi = czero
        eigvecs_full = czero
        eigval_full = zero
        do i=1,nblock
            call csr_to_full(ndim_i, ndim_i, ham_csr(i), ham_full_mpi)
        enddo
        call MPI_ALLREDUCE(ham_full_mpi, ham_full, size(ham_full_mpi), MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        if (myid==master) then
            call full_diag_ham(ndim_i, ham_full, eigval_full, eigvecs_full)
        endif
        call MPI_BARRIER(new_comm, ierror)
        call MPI_BCAST(eigval_full,   size(eigval_full),  MPI_DOUBLE_PRECISION, master, new_comm, ierror)
        call MPI_BCAST(eigvecs_full,  size(eigvecs_full), MPI_DOUBLE_COMPLEX, master, new_comm, ierror)
        eigvals = eigval_full(1:neval)
        eigvecs = eigvecs_full(:,1:nvector)
        deallocate(ham_full)
        deallocate(ham_full_mpi)
        deallocate(eigvecs_full)
        deallocate(eigval_full)
    ! ordinary Lanczos method
    elseif (ed_solver == 1) then
        eigvals = zero
        allocate(eigvecs(nloc,nvector))
        eigvecs = czero
        call diag_ham_lanczos(nblock, end_indx, needed, nloc, neval, maxiter, eigval_tol, ham_csr, eigvals, nvector, eigvecs)
    ! arpack is used
    elseif (ed_solver==2) then
        eigvals = zero
        allocate(eigvecs(nloc,nvector))
        eigvecs = czero
        call diag_ham_arpack(nblock, end_indx, needed, nloc, neval, ham_csr, eigvals, nvector, &
                             eigvecs, ncv, maxiter, eigval_tol, info, actual_step, nconv)
        if (myid==master) then
            print *, "    arpack return info:  ", info
            print *, "    actual arnoldi step:  ", actual_step
            print *, "    number of converged Ritz values:  ", nconv
            print *
        endif
    endif
    call MPI_BARRIER(new_comm, ierror)
    call dealloc_ham_csr(nblock)
    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *
        print *, " fedrixs >>> Done ! Time used: ", time_end-time_begin, "  seconds"
        print *
    endif
    time_begin=time_end

    call write_lowest_eigvals(neval, eigvals)

    if (myid==master) then
        print *, " fedrixs >>> Calculate the density matrix ... "
    endif

    call read_fock_i()
    denmat = czero
    denmat_mpi = czero
    allocate(eigvecs_mpi(ndim_i))
    do k=1, nvector
        eigvecs_mpi = czero
        if (ed_solver==0) then
            eigvecs_mpi = eigvecs(:,k)
        else
            ! send the data
            do i=1,nprocs
               if (myid+1 /= i) then
                   call MPI_ISEND(eigvecs(:,k), nloc, MPI_DOUBLE_COMPLEX, i-1, &
                          i*(10*nprocs)+myid+1, new_comm, request, ierror)
               endif
            enddo
            eigvecs_mpi(end_indx(1,2,myid+1):end_indx(2,2,myid+1)) = eigvecs(:,k)

            ! receive data
            do i=1,nprocs
               if (myid+1 /= i) then
                   call MPI_RECV(eigvecs_mpi(end_indx(1,2,i):end_indx(2,2,i)), &
                        end_indx(2,2,i)-end_indx(1,2,i)+1, MPI_DOUBLE_COMPLEX, &
                       i-1, (myid+1)*(10*nprocs)+i, new_comm, stat, ierror)
               endif
            enddo
        endif
        call MPI_BARRIER(new_comm, ierror)

        do icfg=end_indx(1,1,myid+1), end_indx(2,1,myid+1)
            do j=1,num_val_orbs
                if (btest(fock_i(icfg), j-1)) then
                    denmat_mpi(j,j,k) = denmat_mpi(j,j,k) + conjg(eigvecs_mpi(icfg)) * eigvecs_mpi(icfg)
                endif
            enddo

            if ( abs(eigvecs_mpi(icfg)) < 1E-10 ) cycle
            do i=1,num_val_orbs-1
                do j=i+1,num_val_orbs
                    if (.not. btest(fock_i(icfg), i-1)) cycle
                    old = fock_i(icfg)

                    tot_sign = 1
                    call make_newfock('-', i, old, new, sgn)
                    tot_sign = tot_sign * sgn
                    old = new

                    if (btest(old, j-1)) cycle
                    call make_newfock('+', j,  old, new, sgn)
                    tot_sign = tot_sign * sgn

                    jcfg = binary_search(ndim_i, fock_i, new)

                    if (jcfg == -1) cycle
                    denmat_mpi(i,j,k) = denmat_mpi(i,j,k) + conjg(eigvecs_mpi(icfg)) * eigvecs_mpi(jcfg) * tot_sign
                    denmat_mpi(j,i,k) = denmat_mpi(j,i,k) + conjg(eigvecs_mpi(jcfg)) * eigvecs_mpi(icfg) * tot_sign
                enddo
            enddo
        enddo

        ! dump eigenvectors
        if ( idump ) then
            write(char_I, '(i5)') k
            fname="eigvec."//trim(adjustl(char_I))
            if (myid==master) then
                call write_eigvecs(fname, ndim_i, eigvecs_mpi, eigvals(k))
            endif
        endif
    enddo  ! over k=1, nvector

    call MPI_BARRIER(new_comm, ierror)
    call MPI_ALLREDUCE(denmat_mpi, denmat, size(denmat_mpi), MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
    call write_denmat(nvector, num_val_orbs, denmat)

    call dealloc_fock_i()
    call dealloc_hopping_i()
    call dealloc_coulomb_i()
    deallocate(eigvecs)
    deallocate(eigvecs_mpi)
    deallocate(eigvals)

    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Done ! Time used: ", time_end-time_begin, "  seconds"
        print *
        print *, " fedrixs >>> ED end ! Total time used: ", time_used, "  seconds"
        print *
    endif

    return
end subroutine ed_driver
