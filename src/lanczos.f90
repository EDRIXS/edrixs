!> Lanczos and Krylov-space methods for eigenvalue problems and spectral functions.
!!
!! This module provides three complementary algorithms built on the same
!! three-term Lanczos recurrence  H v_{j+1} = A v_j - alpha_j v_j - beta_{j-1} v_{j-1}:
!!
!!  - diag_ham_lanczos: iterates until the lowest nev eigenvalues converge,
!!    then reconstructs the corresponding eigenvectors by replaying the recurrence.
!!  - build_krylov_mp: runs the recurrence for a fixed number of steps starting
!!    from an externally supplied vector (e.g. T|psi_0>), saving alpha/beta for
!!    the continued-fraction evaluator.
!!  - build_spectrum: evaluates the continued-fraction Green's function
!!    G(omega) = norm / (omega - alpha_1 - beta_1^2/(omega - alpha_2 - ...))
!!    on a user-supplied frequency mesh to yield XAS or RIXS spectral weights.
!!
!! All matrix-vector products are distributed via pspmv_csr so the module
!! scales to any number of MPI ranks.
module m_lanczos

    contains

    !> Parallel Lanczos eigensolver for the lowest nev eigenvalues and eigenvectors.
    !!
    !! Builds a Krylov space using the distributed Hamiltonian ham, starting from
    !! a random initial vector.  Every 10 steps (and at max_niter) the accumulated
    !! tridiagonal matrix is diagonalised with krylov_eigs.  Iteration stops when
    !! the shift in the nev-th eigenvalue falls below eval_tol or the off-diagonal
    !! element beta collapses to zero (exact invariant subspace).
    !!
    !! After convergence the eigenvectors are reconstructed by replaying the
    !! Lanczos recurrence, accumulating the Krylov Ritz vectors.  Only the n_vec
    !! leading eigenvectors are returned (n_vec <= nev).
    !!
    !! @param[in]  nblock     Number of CSR column blocks (equals nprocs)
    !! @param[in]  end_indx   Global index ranges per rank: end_indx(start/end, row/col, rank)
    !! @param[in]  needed     Communication map from get_needed_indx
    !! @param[in]  nloc       Number of local rows / columns
    !! @param[in]  nev        Number of eigenvalues to converge
    !! @param[in]  max_niter  Maximum number of Lanczos iterations
    !! @param[in]  eval_tol   Convergence threshold for eigenvalue shift
    !! @param[in]  ham        Local CSR blocks of the Hamiltonian
    !! @param[out] eval       Converged eigenvalues (length nev), ascending order
    !! @param[in]  n_vec      Number of eigenvectors to return (n_vec <= nev)
    !! @param[out] evec       Local slices of the eigenvectors (nloc x n_vec)
    subroutine diag_ham_lanczos(nblock, end_indx, needed, nloc, nev, max_niter, eval_tol, ham, eval, n_vec, evec)
        use m_constants, only: dp, zero, czero, mystd
        use m_control, only: nprocs , myid, master, new_comm
        use m_types
        use mpi

        implicit none

        integer, intent(in)       :: nblock
        integer, intent(in)       :: end_indx(2,2,nprocs)
        integer, intent(in)       :: needed(nprocs, nprocs)
        integer, intent(in)       :: nloc
        integer, intent(in)       :: nev
        integer, intent(in)       :: max_niter
        real(dp), intent(in)      :: eval_tol
        type (T_csr)              :: ham(nblock)
        real(dp), intent(out)     :: eval(nev)
        integer, intent(in)       :: n_vec
        complex(dp), intent(out)  :: evec(nloc,n_vec)

        real(dp),    allocatable :: temp_vec(:)
        complex(dp), allocatable :: work_vec(:,:)
        real(dp),    allocatable :: alpha(:)
        real(dp),    allocatable :: beta(:)
        real(dp),    allocatable :: krylov_eigval(:)
        real(dp),    allocatable :: krylov_eigvec(:,:)

        integer :: curr
        integer :: prev
        integer :: temp_indx
        integer :: i,j
        integer :: nkrylov
        integer :: ierror

        real(dp) :: res

        complex(dp) :: temp1
        complex(dp) :: temp1_mpi

        allocate(temp_vec(nloc))
        allocate(work_vec(nloc,2))
        allocate(alpha(max_niter))
        allocate(beta(0:max_niter))
        allocate(krylov_eigval(max_niter))
        allocate(krylov_eigvec(max_niter,max_niter))

        temp_vec = zero
        eval     = zero
        work_vec = czero
        alpha    = zero
        beta     = zero

        call random_seed()
        call random_number(temp_vec)
        work_vec(:,1) = temp_vec - 0.5
        if (allocated(temp_vec)) deallocate(temp_vec)

        ! normalize v1
        temp1 = czero
        temp1_mpi = czero
        temp1 = dot_product(work_vec(:,1), work_vec(:,1))
        call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        temp1 = sqrt(real(temp1_mpi))
        work_vec(:,1) = work_vec(:,1) / temp1
        evec(:,1) = work_vec(:,1)

        curr = 1
        prev = 2
        nkrylov = 0

        do i=1,max_niter
            ! matrix-vector product
            call MPI_BARRIER(new_comm, ierror)
            work_vec(:,prev) = -beta(i-1) * work_vec(:,prev)
            call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, work_vec(:,curr), work_vec(:,prev))
            call MPI_BARRIER(new_comm, ierror)
            ! compute alpha_i
            temp1 = czero
            temp1_mpi = czero
            temp1 = dot_product(work_vec(:,curr), work_vec(:,prev))
            call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
            alpha(i) = real(temp1_mpi)

            if (mod(i,10) == 0 .or. i==max_niter) then
                krylov_eigval = zero
                krylov_eigvec = zero
                call krylov_eigs(i, alpha(1:i), beta(1:i-1), krylov_eigval(1:i), krylov_eigvec(1:i,1:i))
                ! check convergence
                if (i>=nev) then
                    res = krylov_eigval(nev) - eval(nev)
                    if (myid==master) then
                        write(mystd,"(a25, i5, E10.2, a5, E10.2)")  "Lanczos iteration:  ", i, abs(res), "-->", eval_tol
                    endif
                    if (abs(res) < eval_tol) then
                        eval(1:nev) = krylov_eigval(1:nev)
                        nkrylov = i
                        EXIT
                    endif
                endif
                eval(1:nev) = krylov_eigval(1:nev)
            endif

            work_vec(:,prev) = work_vec(:,prev) - alpha(i) * work_vec(:,curr)
            ! compute beta_{i+1}
            temp1 = czero
            temp1_mpi = czero
            temp1 = dot_product(work_vec(:,prev), work_vec(:,prev))
            call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
            beta(i) = sqrt(real(temp1_mpi))
            ! beta is too small, exit
            if (abs(beta(i)) < 1E-10) then
                krylov_eigval = zero
                krylov_eigvec = zero
                if (i>=2) then
                    call krylov_eigs(i, alpha(1:i), beta(1:i-1), krylov_eigval(1:i), krylov_eigvec(1:i,1:i))
                    eval(1:nev) = krylov_eigval(1:nev)
                else
                    eval(1) = alpha(1)
                    krylov_eigvec(1,1) = 1.0_dp
                endif
                nkrylov = i
                EXIT
            endif
            work_vec(:,prev) = work_vec(:,prev) / beta(i)
            temp_indx = curr
            curr = prev
            prev = temp_indx
        enddo

        ! build eigenvectors by replaying the Lanczos recurrence
        work_vec = czero
        work_vec(:,1) = evec(:,1)
        evec  = czero
        curr  = 1
        prev  = 2
        do i=1,nkrylov
            do j=1,n_vec
                evec(:,j) = evec(:,j) + krylov_eigvec(i,j) * work_vec(:,curr)
            enddo
            if (i<nkrylov) then
                call MPI_BARRIER(new_comm, ierror)
                work_vec(:,prev) = -alpha(i) * work_vec(:,curr) - beta(i-1) * work_vec(:,prev)
                call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, work_vec(:,curr), work_vec(:,prev))
                call MPI_BARRIER(new_comm, ierror)
                work_vec(:,prev) = work_vec(:,prev) / beta(i)
                temp_indx = curr
                curr = prev
                prev = temp_indx
            endif
        enddo


        if (allocated(work_vec))        deallocate(work_vec)
        if (allocated(alpha))           deallocate(alpha)
        if (allocated(beta))            deallocate(beta)
        if (allocated(krylov_eigval))   deallocate(krylov_eigval)
        if (allocated(krylov_eigvec))   deallocate(krylov_eigvec)

        return
    end subroutine diag_ham_lanczos

    !> Diagonalise a real symmetric tridiagonal matrix using LAPACK dsteqr.
    !!
    !! The tridiagonal matrix is defined by its diagonal alpha and sub-diagonal beta.
    !! On exit eval contains eigenvalues in ascending order and evec the corresponding
    !! orthonormal eigenvectors (columns), which are the Ritz vectors needed to
    !! transform Lanczos basis coordinates back to the original space.
    !!
    !! @param[in]  n     Dimension of the tridiagonal matrix
    !! @param[in]  alpha Diagonal elements (length n)
    !! @param[in]  beta  Sub-diagonal elements (length n-1)
    !! @param[out] eval  Eigenvalues in ascending order (length n)
    !! @param[out] evec  Eigenvectors as columns (n x n)
    subroutine krylov_eigs(n, alpha, beta, eval, evec)
        use m_constants, only: dp, zero
        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: alpha(n)
        real(dp), intent(in) :: beta(n-1)
        real(dp), intent(out) :: eval(n)
        real(dp), intent(out) :: evec(n,n)

        real(dp), allocatable :: tmp_beta(:)
        real(dp), allocatable :: work(:)
        integer :: i
        integer :: info

        allocate(tmp_beta(n-1))
        allocate(work(max(1,2*n-2)))
        eval = alpha
        tmp_beta = beta
        work = zero
        ! init evec to identity matrix
        evec = zero
        do i=1,n
            evec(i,i) = 1.0_dp
        enddo
        call dsteqr('I', n, eval, tmp_beta, evec, n, work, info)
        if (info /= 0 ) then
            print *, "error in dster, info: ", info
            STOP
        endif

        if (allocated(tmp_beta)) deallocate(tmp_beta)
        if (allocated(work))     deallocate(work)

        return
    end subroutine krylov_eigs

    !> Build a Krylov (Lanczos) space starting from an arbitrary initial vector.
    !!
    !! Runs the three-term Lanczos recurrence for up to maxn steps starting from
    !! init_vec.  The actual number of steps taken is returned in neff; it may be
    !! less than maxn if beta collapses (exact Krylov space found).
    !!
    !! The output alpha(1:neff) and beta(0:neff) together with norm define the
    !! continued-fraction representation of the Green's function:
    !!
    !!   G(z) = norm / (z - alpha_1 - beta_1^2/(z - alpha_2 - ...))
    !!
    !! which is evaluated by build_spectrum to obtain XAS or RIXS spectral weights.
    !!
    !! @param[in]  nblock    Number of CSR column blocks
    !! @param[in]  end_indx  Global index ranges per rank
    !! @param[in]  needed    Communication map from get_needed_indx
    !! @param[in]  nloc      Number of local rows / columns
    !! @param[in]  ham       Local CSR blocks of the Hamiltonian
    !! @param[in]  init_vec  Starting vector (length nloc); typically T|psi_0>
    !! @param[in]  maxn      Maximum number of Krylov steps requested
    !! @param[out] neff      Actual number of steps completed
    !! @param[out] alpha     Diagonal elements of the tridiagonal matrix (length maxn)
    !! @param[out] beta      Off-diagonal elements with beta(0)=0 (length 0:maxn)
    !! @param[out] norm      Squared norm of init_vec; serves as the prefactor in G(z)
    subroutine build_krylov_mp(nblock, end_indx, needed, nloc, ham, init_vec, maxn, neff, alpha, beta, norm)
        use m_constants, only: dp, zero, czero, mystd
        use m_control, only: nprocs, myid, master, new_comm
        use m_types
        use mpi

        implicit none

        integer, intent(in)      :: nblock
        integer, intent(in)      :: end_indx(2,2,nprocs)
        integer, intent(in)      :: needed(nprocs, nprocs)
        integer, intent(in)      :: nloc
        type (T_csr)             :: ham(nblock)
        complex(dp), intent(in)  :: init_vec(nloc)
        integer, intent(in)      :: maxn
        integer, intent(out)     :: neff
        real(dp), intent(out)    :: alpha(maxn)
        real(dp), intent(out)    :: beta(0:maxn)
        real(dp), intent(out)    :: norm

        complex(dp), allocatable :: work_vec(:,:)

        integer :: curr
        integer :: prev
        integer :: temp_indx
        integer :: i
        integer :: ierror

        complex(dp) :: temp1
        complex(dp) :: temp1_mpi

        allocate(work_vec(nloc,2))

        work_vec = czero
        alpha    = zero
        beta     = zero
        neff     = 0
        norm     = zero

        work_vec(:,1) = init_vec
        ! normalize V1
        temp1 = czero
        temp1_mpi = czero
        temp1 = dot_product(work_vec(:,1), work_vec(:,1))
        call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        norm = real(temp1_mpi)
        temp1 = sqrt(real(temp1_mpi))
        work_vec(:,1) = work_vec(:,1) / temp1

        curr = 1
        prev = 2

        do i=1,maxn
            if (myid == master .and. mod(i,50)==0) then
                write(mystd,"(a24, i5)")  "Krylov iteration:  ", i
            endif
            neff = i
            ! matrix-vector product
            call MPI_BARRIER(new_comm, ierror)
            work_vec(:,prev) = -beta(i-1) * work_vec(:,prev)
            call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, work_vec(:,curr), work_vec(:,prev))
            call MPI_BARRIER(new_comm, ierror)
            ! compute alpha_i
            temp1 = czero
            temp1_mpi = czero
            temp1 = dot_product(work_vec(:,curr), work_vec(:,prev))
            call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
            alpha(i) = real(temp1_mpi)

            work_vec(:,prev) = work_vec(:,prev) - alpha(i) * work_vec(:,curr)
            ! compute beta_{i+1}
            temp1 = czero
            temp1_mpi = czero
            temp1 = dot_product(work_vec(:,prev), work_vec(:,prev))
            call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
            beta(i) = sqrt(real(temp1_mpi))
            ! beta is too small, exit
            if (abs(beta(i)) < 1E-10) then
                EXIT
            endif
            work_vec(:,prev) = work_vec(:,prev) / beta(i)
            temp_indx = curr
            curr = prev
            prev = temp_indx
        enddo

        if (allocated(work_vec)) deallocate(work_vec)

        return
    end subroutine build_krylov_mp

    !> Evaluate the continued-fraction Green's function on a frequency mesh.
    !!
    !! Given the Lanczos coefficients alpha, beta and the squared norm of the
    !! starting vector, computes the spectral function
    !!
    !!   spec(omega) = -(1/pi) * Im[ norm / (omega + i*gamma + E_gs - CF) ]
    !!
    !! where CF is the continued fraction built bottom-up from level nkrylov to 1.
    !! The ground-state energy e_gs shifts the frequency axis so that the
    !! excitation energy is measured from the ground state.  gamma_final provides
    !! frequency-dependent Lorentzian broadening (final-state lifetime).
    !!
    !! @param[in]  nkrylov      Number of Krylov steps (depth of continued fraction)
    !! @param[in]  alpha        Diagonal Lanczos coefficients (length nkrylov)
    !! @param[in]  beta         Off-diagonal coefficients with beta(0)=0 (length 0:nkrylov)
    !! @param[in]  norm         Squared norm of the starting vector (prefactor)
    !! @param[in]  nw           Number of frequency points
    !! @param[in]  om_mesh      Frequency mesh (eV), length nw
    !! @param[in]  e_gs         Ground-state energy used to shift the frequency axis (eV)
    !! @param[in]  gamma_final  Lorentzian broadening at each frequency point (eV), length nw
    !! @param[out] spec         Spectral function values (length nw)
    subroutine build_spectrum(nkrylov, alpha, beta, norm, nw, om_mesh, e_gs, gamma_final, spec)
        use m_constants, only: dp, czero, pi

        implicit none

        integer, intent(in) :: nkrylov
        real(dp), intent(in) :: alpha(nkrylov)
        real(dp), intent(in) :: beta(0:nkrylov)
        real(dp), intent(in) :: norm
        integer,  intent(in) :: nw
        real(dp), intent(in) :: om_mesh(nw)
        real(dp), intent(in) :: e_gs
        real(dp), intent(in) :: gamma_final(nw)
        real(dp), intent(out) :: spec(nw)

        integer :: i
        complex(dp) :: tmp_vec(nw)

        tmp_vec = czero
        do i=nkrylov, 2, -1
            tmp_vec = beta(i-1)**2 / (om_mesh + dcmplx(0.0, 1.0)*gamma_final + e_gs - alpha(i) - tmp_vec)
        enddo
        tmp_vec = 1.0/(om_mesh + dcmplx(0.0,1.0) * gamma_final + e_gs - alpha(1) - tmp_vec)

        spec = -1.0/pi * aimag(tmp_vec) * norm

        return
    end subroutine build_spectrum

end module m_lanczos
