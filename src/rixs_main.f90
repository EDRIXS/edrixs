program rixs_main
    use m_control, only: master, origin_myid, origin_nprocs, origin_comm
    use m_control, only: myid, nprocs, new_comm, ndim_n, ndim_i, ndim_f, ndim_n_nocore, num_core_orbs
    use m_global, only: dealloc_fock_i, dealloc_fock_n, dealloc_fock_f
    use mpi 

    implicit none
 
    integer :: ierror
    integer :: color
    integer :: key
    integer :: min_dim

    origin_comm = MPI_COMM_WORLD
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(origin_comm, origin_myid,   ierror)
    call MPI_COMM_SIZE(origin_comm, origin_nprocs, ierror)
    call MPI_BARRIER(origin_comm, ierror)

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
    call MPI_FINALIZE(ierror)

end program rixs_main
