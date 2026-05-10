!> Full diagonalisation of a complex Hermitian matrix via LAPACK ZHEEV.
!!
!! Used when the Hilbert-space dimension is small enough (ndim_i < min_ndim)
!! that storing and diagonalising the dense matrix is cheaper than running
!! the iterative Lanczos solver.  All eigenvalues and eigenvectors are
!! returned; the caller selects the neval lowest ones.
!!
!! @param[in]  ndim  Matrix dimension
!! @param[in]  ham   Complex Hermitian matrix to diagonalise (ndim x ndim)
!! @param[out] eval  Eigenvalues in ascending order (length ndim)
!! @param[out] evec  Corresponding eigenvectors as columns (ndim x ndim)
subroutine full_diag_ham(ndim, ham, eval, evec)
    use m_constants, only: dp, zero

    implicit none

    integer,     intent(in)  :: ndim
    complex(dp), intent(in)  :: ham(ndim, ndim)
    real(dp),    intent(out) :: eval(ndim)
    complex(dp), intent(out) :: evec(ndim, ndim)

    integer :: info
    integer :: lwork
    integer :: lrwork
    real(dp),    allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)

    lwork  = 2 * ndim - 1
    lrwork = 3 * ndim - 2

    allocate(work(lwork))
    allocate(rwork(lrwork))

    eval = zero
    evec = ham

    call ZHEEV('V', 'U', ndim, evec, ndim, eval, work, lwork, rwork, info)

    if (allocated(work))  deallocate(work)
    if (allocated(rwork)) deallocate(rwork)

    return
end subroutine full_diag_ham
