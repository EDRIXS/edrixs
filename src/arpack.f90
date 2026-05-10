!> Parallel ARPACK eigensolver for the lowest nev eigenvalues and eigenvectors.
!!
!! Uses the parallel-ARPACK reverse-communication routine pznaupd to build an
!! Implicitly Restarted Arnoldi (IRAM) factorisation of the distributed sparse
!! Hamiltonian.  Each time pznaupd requests a matrix-vector product (ido=-1 or
!! ido=1), we apply pspmv_csr to the input slice in workd(ipntr(1):) and write
!! the result into workd(ipntr(2):).  When pznaupd signals convergence
!! (any other value of ido), we call pzneupd to extract the eigenpairs.
!!
!! The default settings request the nev eigenvalues with smallest real part
!! (which='SR') in regular mode (mode=1, no shift-invert), which is appropriate
!! for finding the ground state of a Hermitian Hamiltonian.  Although the
!! Hamiltonian is Hermitian, the complex (znaupd/zneupd) interface is used
!! because the Hamiltonian elements are stored as double-complex.
!!
!! ARPACK requires ncv > nev + 1 working Arnoldi vectors; ncv ~ 2*nev or larger
!! gives faster convergence but uses more memory.  This is the recommended
!! eigensolver (ed_solver=2) for SIAM and large multi-shell ED problems.
!!
!! @param[in]  nblock      Number of CSR column blocks (equals nprocs)
!! @param[in]  end_indx    Global index ranges per rank: end_indx(start/end, row/col, rank)
!! @param[in]  needed      Communication map from get_needed_indx
!! @param[in]  nloc        Number of local rows / columns owned by this rank
!! @param[in]  nev         Number of eigenvalues to compute
!! @param[in]  ham         Local CSR blocks of the Hamiltonian
!! @param[out] eval        Converged eigenvalues, ascending order (length nev)
!! @param[in]  n_vecs      Number of eigenvectors to return (n_vecs <= nev)
!! @param[out] evec        Local slices of the eigenvectors (nloc x n_vecs)
!! @param[in]  ncv         Size of the Arnoldi basis (ncv > nev + 1; ncv ~ 2*nev recommended)
!! @param[in]  maxiter     Maximum number of IRAM iterations
!! @param[in]  tol         Relative convergence tolerance for Ritz values
!! @param[out] info        ARPACK return code: 0 = converged, < 0 = error, > 0 = warning
!! @param[out] actual_step Actual number of IRAM iterations performed (iparam(3))
!! @param[out] nconv       Number of converged Ritz values (iparam(5))
subroutine diag_ham_arpack(nblock, end_indx, needed, nloc, nev, ham, eval, &
                   n_vecs, evec, ncv, maxiter, tol, info, actual_step, nconv)

    use m_constants, only: dp, czero
    use m_control,   only: nprocs, new_comm
    use m_types
    use mpi

    implicit none

    integer, intent(in)       :: nblock
    integer, intent(in)       :: end_indx(2,2,nprocs)
    integer, intent(in)       :: needed(nprocs, nprocs)
    integer, intent(in)       :: nloc
    integer, intent(in)       :: nev
    type (T_csr)              :: ham(nblock)
    real(dp), intent(out)     :: eval(nev)
    integer, intent(in)       :: n_vecs
    complex(dp), intent(out)  :: evec(nloc,n_vecs)
    integer, intent(in)       :: ncv
    integer, intent(in)       :: maxiter
    real(dp), intent(in)      :: tol
    integer, intent(out)      :: info
    integer, intent(out)      :: actual_step
    integer, intent(out)      :: nconv

    integer      :: ido          ! Reverse-communication flag (set by pznaupd)
    integer      :: ishfts        ! 1 = exact shifts (recommended)
    integer      :: mode          ! 1 = regular (A x = lambda x, no shift-invert)
    integer      :: lworkl        ! Size of workl, formula from ARPACK docs
    integer      :: iparam(11)    ! ARPACK parameter array (see pznaupd docs)
    integer      :: ipntr(14)     ! ARPACK pointer array into workd / workl
    integer      :: sorted_indx(nev)
    integer      :: i

    logical      :: rvec          ! .true. = compute eigenvectors as well as eigenvalues
    character(1) :: bmat          ! 'I' = standard (not generalised) eigenproblem
    character(2) :: which         ! 'SR' = smallest real part (lowest energy)
    complex(dp)  :: sigma         ! Shift parameter (unused in mode=1)

    real(dp)   , allocatable  :: rwork(:)
    logical    , allocatable  :: selec(:)
    complex(dp), allocatable  :: resid(:)    ! Residual vector for restart
    complex(dp), allocatable  :: v(:,:)      ! Arnoldi basis vectors
    complex(dp), allocatable  :: workd(:)    ! Reverse-communication work space
    complex(dp), allocatable  :: workl(:)    ! Internal ARPACK work space
    complex(dp), allocatable  :: workev(:)
    complex(dp), allocatable  :: d(:)        ! Computed Ritz values

    lworkl  = 3*ncv**2 + 5*ncv

    allocate(rwork(ncv))
    allocate(selec(ncv))
    allocate(resid(nloc))
    allocate(v(nloc,ncv))
    allocate(workd(3*nloc))
    allocate(workl(lworkl))
    allocate(workev(3*ncv))
    allocate(d(ncv))

    ishfts    = 1
    mode      = 1
    iparam(1) = ishfts
    iparam(3) = maxiter
    iparam(7) = mode
    bmat      = 'I'
    which     = 'SR'
    ido       = 0
    info      = 0
    selec     = .true.

    ! Reverse-communication loop: pznaupd asks us to apply H to a vector
    ! whenever ido = -1 or 1; any other value means "done, call pzneupd".
    do while (.true.)
        call pznaupd(new_comm, ido, bmat, nloc, which, nev, tol, resid, ncv, &
                     v, nloc, iparam, ipntr, workd, workl, lworkl, rwork, info)

        if (ido .eq. -1 .or. ido .eq. 1) then
            workd(ipntr(2):ipntr(2)+nloc-1) = czero
            call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, &
                                                     workd(ipntr(1)), workd(ipntr(2)))
        else
            EXIT
        endif
    enddo
    actual_step = iparam(3)
    nconv       = iparam(5)

    if (info .lt. 0) then
        return
    else
        ! Extract eigenvalues and eigenvectors from the converged Arnoldi factorisation
        rvec = .true.
        call pzneupd(new_comm, rvec, 'A', selec, d, v, nloc, sigma, workev, &
                          bmat, nloc, which, nev, tol, resid, ncv, v, nloc, iparam, &
                          ipntr, workd, workl, lworkl, rwork, info)
        if (info .ne. 0) then
            return
        endif
        ! ARPACK does not guarantee any ordering of Ritz values; sort ascending
        eval = real(d(1:nev))
        do i=1,nev
            sorted_indx(i) = i
        enddo
        call sort_eigvals(nev, eval, sorted_indx)
        do i=1,n_vecs
            evec(:,i) = v(:,sorted_indx(i))
        enddo
    endif

    return
end subroutine diag_ham_arpack
