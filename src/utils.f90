!> Binary search for value m in a sorted integer array vec of length n.
!!
!! Returns the 1-based index of m in vec, or -1 if m is not found.
!! The array must be sorted in ascending order.
!!
!! @param[in] n    Length of vec
!! @param[in] vec  Sorted array of 64-bit integers
!! @param[in] m    Value to search for
!! @return         1-based index of m in vec, or -1 if not found
function binary_search(n, vec, m) result(val)
    use m_constants, only: dp

    implicit none

    integer,    intent(in) :: n
    integer(dp), intent(in) :: vec(n)
    integer(dp), intent(in) :: m

    integer :: val

    integer :: left, right, middle

    val   = -1
    left  = 1
    right = n
    do while (left <= right)
        middle = left + (right - left) / 2
        if (vec(middle) < m) then
            left = middle + 1
        elseif (vec(middle) > m) then
            right = middle - 1
        else
            val = middle
            EXIT
        endif
    enddo

    return
end function binary_search

!> Sort eigenvalues into ascending order, keeping the index permutation.
!!
!! Uses insertion sort, which is efficient for nearly-sorted small arrays
!! as typically encountered after Lanczos convergence.
!!
!! @param[in]    n           Number of eigenvalues
!! @param[inout] eigvals     On entry: unsorted eigenvalues.  On exit: sorted ascending.
!! @param[inout] sorted_indx On entry: initial index array.  On exit: permutation that sorts eigvals.
subroutine sort_eigvals(n, eigvals, sorted_indx)
    use m_constants, only: dp, zero
    implicit none

    integer,  intent(in)    :: n
    real(dp), intent(inout) :: eigvals(n)
    integer,  intent(inout) :: sorted_indx(n)

    real(dp) :: tmp_eigval(n)
    integer  :: tmp_indx(n)
    integer  :: i, j

    tmp_eigval = eigvals
    tmp_indx   = sorted_indx

    eigvals     = zero
    sorted_indx = 0
    eigvals(1)     = tmp_eigval(1)
    sorted_indx(1) = tmp_indx(1)
    do i = 2, n
        j = i - 1
        do while (tmp_eigval(i) < eigvals(j))
            eigvals(j+1)     = eigvals(j)
            sorted_indx(j+1) = sorted_indx(j)
            j = j - 1
            if (j == 0) EXIT
        enddo
        eigvals(j+1)     = tmp_eigval(i)
        sorted_indx(j+1) = tmp_indx(i)
    enddo

    return
end subroutine sort_eigvals

!> Partition m rows and n columns evenly across nprocs MPI ranks.
!!
!! Returns a 3-D index array end_indx(2,2,nprocs) where
!!   end_indx(1,1,p) / end_indx(2,1,p) are the first/last global row indices for rank p,
!!   end_indx(1,2,p) / end_indx(2,2,p) are the first/last global column indices for rank p.
!! The last rank absorbs any remainder from integer division.
!!
!! @param[in]  nprocs    Total number of MPI ranks
!! @param[in]  m         Total number of rows to partition
!! @param[in]  n         Total number of columns to partition
!! @param[out] end_indx  Index ranges: end_indx(start/end, row/col, rank)
subroutine partition_task(nprocs, m, n, end_indx)
    implicit none

    integer, intent(in)  :: nprocs
    integer, intent(in)  :: m
    integer, intent(in)  :: n
    integer, intent(out) :: end_indx(2, 2, nprocs)

    integer :: step, i

    end_indx = 0

    step = m / nprocs
    do i = 1, nprocs-1
        end_indx(1,1,i) = (i-1)*step + 1
        end_indx(2,1,i) =     i*step
    enddo
    end_indx(1,1,nprocs) = (nprocs-1)*step + 1
    end_indx(2,1,nprocs) = m

    step = n / nprocs
    do i = 1, nprocs-1
        end_indx(1,2,i) = (i-1)*step + 1
        end_indx(2,2,i) =     i*step
    enddo
    end_indx(1,2,nprocs) = (nprocs-1)*step + 1
    end_indx(2,2,nprocs) = n

    return
end subroutine partition_task

!> Build the global communication map for a distributed sparse matrix.
!!
!! For each pair (myid, i), sets needed(myid+1, i) = 1 if this rank's local
!! Hamiltonian has a non-zero off-diagonal block targeting column-block i.
!! The map is then reduced with MPI_ALLREDUCE so every rank knows which
!! column vectors it must send to others and which it will receive.
!!
!! @param[in]  nblock    Number of CSR column blocks
!! @param[in]  ham_csr   Array of local CSR blocks
!! @param[out] needed    Communication map: needed(i,j)=1 means rank i needs data from rank j
subroutine get_needed_indx(nblock, ham_csr, needed)
    use m_control, only: myid, nprocs, new_comm
    use m_types
    use mpi
    implicit none

    integer, intent(in)  :: nblock
    type(T_csr)          :: ham_csr(nblock)
    integer, intent(out) :: needed(nprocs, nprocs)

    integer :: i, ierror
    integer :: needed_mpi(nprocs, nprocs)

    needed     = 0
    needed_mpi = 0
    if (nblock == nprocs) then
        do i = 1, nprocs
            if (myid+1 /= i .and. ham_csr(i)%nnz > 0) then
                needed(myid+1, i) = 1
            endif
        enddo
        call MPI_BARRIER(new_comm, ierror)
        call MPI_ALLREDUCE(needed, needed_mpi, size(needed), MPI_INTEGER, MPI_SUM, new_comm, ierror)
        call MPI_BARRIER(new_comm, ierror)
        needed = needed_mpi
    endif

    return
end subroutine get_needed_indx

!> Sum the total number of non-zero elements across all local CSR blocks.
!!
!! Each rank sums its own blocks, then MPI_ALLREDUCE gives the global total.
!!
!! @param[in]  nblock        Number of local CSR blocks
!! @param[in]  ham_csr       Array of local CSR blocks
!! @param[out] num_nonzeros  Global total number of non-zero elements
subroutine get_number_nonzeros(nblock, ham_csr, num_nonzeros)
    use m_constants, only: dp
    use m_control,   only: new_comm
    use m_types
    use mpi
    implicit none

    integer,     intent(in)  :: nblock
    type(T_csr)              :: ham_csr(nblock)
    integer(dp), intent(out) :: num_nonzeros

    integer :: i, ierror
    integer(dp) :: num_nonzeros_mpi

    num_nonzeros     = 0_dp
    num_nonzeros_mpi = 0_dp
    do i = 1, nblock
        num_nonzeros = num_nonzeros + ham_csr(i)%nnz
    enddo
    call MPI_BARRIER(new_comm, ierror)
    call MPI_ALLREDUCE(num_nonzeros, num_nonzeros_mpi, 1, MPI_INTEGER8, MPI_SUM, new_comm, ierror)
    call MPI_BARRIER(new_comm, ierror)
    num_nonzeros = num_nonzeros_mpi

    return
end subroutine get_number_nonzeros

!> Compute Boltzmann thermal occupation probabilities for a set of eigenstates.
!!
!! Shifts energies by their minimum before exponentiation to avoid overflow,
!! then normalises by the partition function.
!!
!! @param[in]  n       Number of eigenstates
!! @param[in]  eigval  Eigenvalues (eV)
!! @param[in]  beta    Inverse temperature beta = 1/(k_B T) in 1/eV
!! @param[out] prob    Boltzmann weights, normalised so sum(prob) = 1
subroutine get_prob(n, eigval, beta, prob)
    use m_constants
    implicit none

    integer,  intent(in)  :: n
    real(dp), intent(in)  :: eigval(n)
    real(dp), intent(in)  :: beta
    real(dp), intent(out) :: prob(n)

    real(dp) :: rtemp(n)
    real(dp) :: norm
    integer  :: i

    rtemp = eigval - minval(eigval)
    norm  = 0.0
    do i = 1, n
        norm = norm + exp(-beta * rtemp(i))
    enddo

    prob = exp(-beta * rtemp) / norm

    return
end subroutine get_prob
