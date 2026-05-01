!> Python API entry points exposing the four EDRIXS solvers via f2py.
!!
!! Each *_fsolver subroutine accepts an MPI communicator, rank, and size from
!! the calling Python layer (typically obtained from mpi4py) and handles the
!! rank-reduction logic: if the user launched more MPI processes than there are
!! Hilbert-space basis states, the excess ranks are split off into a second
!! communicator with MPI_COMM_SPLIT and excluded from the actual computation.
!! This prevents zero-row partitions that would crash the CSR builder.
!!
!! After rank reduction each active rank calls the corresponding *_driver
!! subroutine, then all ranks (including idle ones) synchronise on the original
!! communicator via MPI_BARRIER before returning.

!> Python-callable wrapper for the ED solver.
!!
!! Accepts the MPI communicator from the Python side, performs rank reduction
!! if nprocs > ndim_i, then calls ed_driver on active ranks.
!!
!! @param[in] comm      MPI communicator handle (from mpi4py)
!! @param[in] my_id     Rank of the calling process in comm
!! @param[in] num_procs Total number of processes in comm
subroutine ed_fsolver(comm, my_id, num_procs)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_i
    use m_global, only: dealloc_fock_i
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs

    integer :: ierror
    integer :: color
    integer :: key

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config()
    ! read fock to know the dimension of the Hamiltonian
    call read_fock_i()
    call dealloc_fock_i()
    if (ndim_i < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than ndim_i: ", ndim_i
            print *, " fedrixs >>> Only ", ndim_i, " processors will really work!"
        endif
        if (origin_myid < ndim_i) then
            color = 1
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < ndim_i) then
        call ed_driver()
    endif

    call MPI_BARRIER(origin_comm, ierror)
    return
end subroutine ed_fsolver

!> Python-callable wrapper for the XAS solver.
!!
!! Performs rank reduction if nprocs > min(ndim_i, ndim_n), then calls
!! xas_driver on active ranks.
!!
!! @param[in] comm      MPI communicator handle (from mpi4py)
!! @param[in] my_id     Rank of the calling process in comm
!! @param[in] num_procs Total number of processes in comm
subroutine xas_fsolver(comm, my_id, num_procs)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n, ndim_i, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_i, dealloc_fock_n
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs

    integer :: ierror
    integer :: color
    integer :: key
    integer :: min_dim

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config()
    call read_fock_i()
    call dealloc_fock_i()
    call read_fock_n()
    call dealloc_fock_n()
    ndim_n = ndim_n_nocore * num_core_orbs
    min_dim = min(ndim_i, ndim_n)
    if (min_dim < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than min(ndim_i, ndim_n): ", ndim_i, ndim_n
            print *, " fedrixs >>> Only ", min_dim, " processors will really work!"
        endif
        if (origin_myid < min_dim) then
            color = 1
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < min_dim) then
       call xas_driver()
    endif

    call MPI_BARRIER(origin_comm, ierror)

    return
end subroutine xas_fsolver

!> Python-callable wrapper for the RIXS solver.
!!
!! Performs rank reduction if nprocs > min(ndim_i, ndim_n, ndim_f), then calls
!! rixs_driver on active ranks.
!!
!! @param[in] comm      MPI communicator handle (from mpi4py)
!! @param[in] my_id     Rank of the calling process in comm
!! @param[in] num_procs Total number of processes in comm
subroutine rixs_fsolver(comm, my_id, num_procs)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n, ndim_i, ndim_f, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_i, dealloc_fock_n, dealloc_fock_f
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs

    integer :: ierror
    integer :: color
    integer :: key
    integer :: min_dim

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs

    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config()
    call read_fock_i()
    call dealloc_fock_i()
    call read_fock_n()
    call dealloc_fock_n()
    call read_fock_f()
    call dealloc_fock_f()
    ndim_n = ndim_n_nocore * num_core_orbs
    min_dim = min(ndim_i, ndim_n, ndim_f)
    if (min_dim < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than min(ndim_i, ndim_n, ndim_f): ", ndim_i, ndim_n, ndim_f
            print *, " fedrixs >>> Only ", min_dim, " processors will really work!"
        endif
        if (origin_myid < min_dim) then
            color = 1
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif

    if (origin_myid < min_dim) then
        call rixs_driver()
    endif

    call MPI_BARRIER(origin_comm, ierror)
    return
end subroutine rixs_fsolver

!> Python-callable wrapper for the operator-average solver.
!!
!! Performs rank reduction if nprocs > ndim_i, then calls opavg_driver on
!! active ranks.
!!
!! @param[in] comm      MPI communicator handle (from mpi4py)
!! @param[in] my_id     Rank of the calling process in comm
!! @param[in] num_procs Total number of processes in comm
subroutine opavg_fsolver(comm, my_id, num_procs)
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_i
    use m_global, only: dealloc_fock_i
    use mpi

    implicit none

    integer, intent(in) :: comm
    integer, intent(in) :: my_id
    integer, intent(in) :: num_procs

    integer :: ierror
    integer :: color
    integer :: key

!F2PY intent(in) comm
!F2PY intent(in) my_id
!F2PY intent(in) num_procs
    origin_comm = comm
    origin_myid = my_id
    origin_nprocs = num_procs

    call config()
    call read_fock_i()
    call dealloc_fock_i()
    if (ndim_i < origin_nprocs) then
        if (origin_myid==master) then
            print *, " fedrixs >>> Warning: number of CPU processors ", origin_nprocs, &
                     "is larger than ndim_i: ", ndim_i
            print *, " fedrixs >>> Only ", ndim_i, " processors will really work!"
        endif
        if (origin_myid < ndim_i) then
            color = 1
            key = origin_myid
        else
            color = 2
            key = origin_myid
        endif
        call MPI_COMM_SPLIT(origin_comm, color, key, new_comm, ierror)
        call MPI_COMM_RANK(new_comm, myid, ierror)
        call MPI_COMM_SIZE(new_comm, nprocs, ierror)
    else
        myid = origin_myid
        nprocs = origin_nprocs
        new_comm = origin_comm
    endif
    print *, " fedrixs >>> ", origin_myid, origin_nprocs, myid, nprocs

    if (origin_myid < ndim_i) then
        call opavg_driver()
    endif

    call MPI_BARRIER(origin_comm, ierror)

    return
end subroutine opavg_fsolver
