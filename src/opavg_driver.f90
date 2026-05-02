!> Operator-average driver: compute <psi_k|H|psi_k> for stored eigenstates.
!!
!! Reads the initial-state Hamiltonian and applies it to each eigenvector
!! stored on disk by ed.x.  For state k the expectation value is
!!
!!   opavg(k) = <psi_k| H |psi_k>
!!
!! computed as a distributed dot product v^H (H v) using pspmv_csr.
!! Results are written to opavg.dat via write_opavg.
!!
!! This is useful for verifying the ED eigenvalues and for computing
!! operator averages when the operator shares the same Fock basis as H_i.
subroutine opavg_driver()
    use m_constants
    use m_control
    use m_global
    use mpi

    implicit none

    integer :: mblock
    integer :: nblock
    integer :: nloc
    integer :: needed(nprocs,nprocs)
    integer :: end_indx(2,2,nprocs)
    integer :: ierror
    integer :: igs
    integer(dp) :: num_of_nonzeros

    real(dp) :: rtemp
    real(dp) :: time_begin
    real(dp) :: time_end
    real(dp) :: time_used
    real(dp) :: eigvals(num_gs)

    complex(dp) :: omega
    complex(dp) :: temp1
    complex(dp) :: temp1_mpi
    complex(dp) :: opavg(num_gs)
    complex(dp), allocatable :: v_vector(:)
    complex(dp), allocatable :: w_vector(:)
    complex(dp), allocatable :: v_vector_mpi(:)

    character(len=20) :: fname
    character(len=10) :: char_I

    time_begin = 0.0_dp
    time_end   = 0.0_dp
    time_used  = 0.0_dp
    call cpu_time(time_begin)

    if (myid == master) then
        print *, "--------------------------------------------"
        print *, " fedrixs >>> OPAVG begin ... "
        print *
    endif
    call read_hopping_i()
    call read_coulomb_i()
    call read_fock_i()

    if (myid == master) then
        write(mystd,"(a20, i15)")   "num_val_orbs:  ", num_val_orbs
        write(mystd,"(a20, i15)")   "nhopp_i:       ", nhopp_i
        write(mystd,"(a20, i15)")   "ncoul_i:       ", ncoul_i
        write(mystd,"(a20, i15)")   "ndim_i:        ", ndim_i
        write(mystd,"(a20, i15)")   "nprocs:        ", nprocs
        write(mystd,"(a20, i15)")   "num_gs:        ", num_gs
        print *
    endif

    if (myid==master) then
        print *, " fedrixs >>> Building operator  ..."
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
    call dealloc_fock_i()
    call get_needed_indx(nblock, ham_csr, needed)
    call get_number_nonzeros(nblock, ham_csr, num_of_nonzeros)
    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Number of nonzero elements of the Hamiltonian", num_of_nonzeros
        print *, " fedrixs >>> Done !"
        print *, " fedrixs >>> Time used: ", time_end - time_begin, "  seconds"
        print *
    endif

    eigvals=zero
    opavg=czero
    allocate(v_vector(nloc))
    allocate(w_vector(nloc))
    time_begin = time_end
    do igs=1, num_gs
        if (myid == master) then
            print *, " fedrixs >>> For state: ", igs
        endif
        allocate(v_vector_mpi(ndim_i))
        v_vector = czero
        w_vector = czero
        v_vector_mpi = czero
        write(char_I, '(i5)') igs
        fname="eigvec."//trim(adjustl(char_I))
        call read_eigvecs(fname, ndim_i, v_vector_mpi, eigvals(igs))
        v_vector = v_vector_mpi(end_indx(1,2,myid+1): end_indx(2,2,myid+1))
        deallocate(v_vector_mpi)
        call MPI_BARRIER(new_comm, ierror)
        call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham_csr, v_vector, w_vector)
        call MPI_BARRIER(new_comm, ierror)
        temp1 = czero
        temp1_mpi = czero
        temp1 = dot_product(v_vector, w_vector)
        call MPI_ALLREDUCE(temp1, temp1_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        opavg(igs) = temp1_mpi
    enddo

    deallocate(v_vector)
    deallocate(w_vector)

    call write_opavg(num_gs, eigvals, opavg)

    call cpu_time(time_end)
    time_used = time_used + time_end - time_begin
    if (myid==master) then
        print *, " fedrixs >>> Done ! Time used: ", time_end-time_begin, "  seconds"
        print *
        print *, " fedrixs >>> OPAVG end ! Total time used: ", time_used, "  seconds"
        print *
    endif

    return
end subroutine opavg_driver
