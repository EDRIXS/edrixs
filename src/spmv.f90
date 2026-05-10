!> Parallel distributed sparse matrix–vector product using MPI.
!!
!! The global matrix is stored as nblock CSR blocks distributed across MPI
!! ranks.  Each rank owns a contiguous row-block and one column-block (the
!! diagonal block).  Off-diagonal column-blocks reside on other ranks and
!! must be communicated before the corresponding local multiply can proceed.
!!
!! Communication pattern:
!!  1. Each rank non-blocking sends its local vector slice to every rank that
!!     needs it (as indicated by the needed(:,:) map).
!!  2. While sends are in flight the diagonal block multiply is done locally.
!!  3. Off-diagonal block vectors are received one by one and each multiply
!!     is accumulated into vec_out immediately after receipt.
!!  4. MPI_WAITALL ensures all sends have completed before returning.
!!
!! @param[in]    comm      MPI communicator
!! @param[in]    nblock    Number of CSR column-blocks (equals nprocs in the distributed case)
!! @param[in]    end_indx  Global index ranges for each rank: end_indx(start/end, row/col, rank)
!! @param[in]    needed    needed(i,j)=1 means rank i needs the column vector from rank j
!! @param[in]    mloc      Number of local rows
!! @param[in]    nloc      Number of local columns
!! @param[inout] ham       Array of nblock CSR blocks for this rank's rows
!! @param[in]    vec_in    Local input vector (length nloc)
!! @param[inout] vec_out   Local output vector (length mloc); contributions are accumulated
subroutine pspmv_csr(comm, nblock, end_indx, needed, mloc, nloc, ham, vec_in, vec_out)
    use m_constants, only: dp, czero
    use m_control,   only: myid, nprocs
    use m_types,     only: T_csr
    use mpi

    implicit none

    integer, intent(in)        :: comm
    integer, intent(in)        :: nblock
    integer, intent(in)        :: end_indx(2, 2, nprocs)
    integer, intent(in)        :: needed(nprocs, nprocs)
    integer, intent(in)        :: mloc
    integer, intent(in)        :: nloc
    type(T_csr)                :: ham(nblock)
    complex(dp), intent(in)    :: vec_in(nloc)
    complex(dp), intent(inout) :: vec_out(mloc)

    integer :: i, ireq, ierror
    integer :: request(nprocs)
    integer :: send_stat(MPI_STATUS_SIZE, nprocs)
    integer :: stat(MPI_STATUS_SIZE)
    integer :: other_size, max_size
    complex(dp), allocatable :: tmp_vec(:)

    ireq = 0
    do i = 1, nprocs
        if (needed(i, myid+1) == 1) then
            ireq = ireq + 1
            call MPI_ISEND(vec_in, nloc, MPI_DOUBLE_COMPLEX, i-1, &
                           i*(10*nprocs)+myid+1, comm, request(ireq), ierror)
        endif
    enddo

    ! Diagonal block: no communication needed
    call matvec_csr(mloc, nloc, ham(myid+1), vec_in, vec_out)

    ! Off-diagonal blocks: receive remote column slices and multiply
    if (nblock == nprocs .and. nprocs > 1) then
        max_size = end_indx(2,2,1) - end_indx(1,2,1) + 1
        do i = 2, nprocs
            if ((end_indx(2,2,i) - end_indx(1,2,i) + 1) > max_size) &
                max_size = end_indx(2,2,i) - end_indx(1,2,i) + 1
        enddo
        allocate(tmp_vec(max_size))
        do i = 1, nprocs
            if (needed(myid+1, i) == 1) then
                other_size = end_indx(2,2,i) - end_indx(1,2,i) + 1
                tmp_vec = czero
                call MPI_RECV(tmp_vec(1:other_size), other_size, MPI_DOUBLE_COMPLEX, &
                              i-1, (myid+1)*(10*nprocs)+i, comm, stat, ierror)
                call matvec_csr(mloc, other_size, ham(i), tmp_vec(1:other_size), vec_out)
            endif
        enddo
        if (allocated(tmp_vec)) deallocate(tmp_vec)
    endif

    call MPI_WAITALL(ireq, request, send_stat, ierror)

    return
end subroutine pspmv_csr

!> Local (single-rank) sparse matrix–vector product for one CSR block.
!!
!! Computes vec_out += mat * vec_in where mat is stored in CSR format.
!! The col_shift offset in mat maps local column indices back to the global
!! numbering so that vec_in can be indexed correctly.
!!
!! @param[in]    m        Number of rows of mat
!! @param[in]    n        Length of vec_in
!! @param[in]    mat      CSR matrix block
!! @param[in]    vec_in   Input vector (length n)
!! @param[inout] vec_out  Output vector (length m); result is accumulated into this
subroutine matvec_csr(m, n, mat, vec_in, vec_out)
    use m_constants, only: dp
    use m_types

    implicit none

    integer,     intent(in)    :: m
    integer,     intent(in)    :: n
    type(T_csr), intent(in)    :: mat
    complex(dp), intent(in)    :: vec_in(n)
    complex(dp), intent(inout) :: vec_out(m)

    complex(dp) :: tmp
    integer :: i, j, begin_col, end_col

    do i = 1, mat%m
        begin_col = mat%iaa(i)
        end_col   = mat%iaa(i+1) - 1
        if (begin_col > end_col) cycle
        tmp = 0
        do j = begin_col, end_col
            tmp = tmp + mat%aa(j) * vec_in(mat%jaa(j) - mat%col_shift)
        enddo
        vec_out(i) = vec_out(i) + tmp
    enddo

    return
end subroutine matvec_csr
