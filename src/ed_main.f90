program ed_main
    use m_control, only: origin_myid, origin_nprocs, origin_comm
    use m_control, only: master, myid, nprocs, new_comm, ndim_i
    use m_global, only: dealloc_fock_i
    use mpi 

    implicit none
 
    integer :: ierror
    integer :: color
    integer :: key

    origin_comm = MPI_COMM_WORLD
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(origin_comm, origin_myid,   ierror)
    call MPI_COMM_SIZE(origin_comm, origin_nprocs, ierror)
    call MPI_BARRIER(origin_comm, ierror)

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
    call MPI_FINALIZE(ierror)

end program ed_main
