!> Global allocatable arrays for Hamiltonian data, Fock bases, and Krylov vectors.
!!
!! All arrays are module-level with the save attribute. They are allocated and
!! deallocated through the paired alloc_* / dealloc_* subroutines defined here,
!! which read the required sizes from m_control.
module m_global
    use m_constants, only: dp
    use m_types

    implicit none

    !> Non-zero hopping elements of the initial (no core-hole) Hamiltonian.
    type(T_tensor2), public, allocatable, save :: hopping_i(:)
    !> Non-zero Coulomb elements of the initial Hamiltonian.
    type(T_tensor4), public, allocatable, save :: coulomb_i(:)
    !> Non-zero hopping elements of the intermediate (core-hole) Hamiltonian.
    type(T_tensor2), public, allocatable, save :: hopping_n(:)
    !> Non-zero Coulomb elements of the intermediate Hamiltonian.
    type(T_tensor4), public, allocatable, save :: coulomb_n(:)

    !> Non-zero elements of the XAS photon-absorption transition operator.
    type(T_tensor2), public, allocatable, save :: transop_xas(:)
    !> Non-zero elements of the RIXS absorption transition operator (initial -> intermediate).
    type(T_tensor2), public, allocatable, save :: transop_rixs_i(:)
    !> Non-zero elements of the RIXS emission transition operator (intermediate -> final).
    type(T_tensor2), public, allocatable, save :: transop_rixs_f(:)

    !> Fock basis for the initial Hilbert space (length ndim_i).
    !! Each element is an integer whose bits encode orbital occupations.
    integer(dp), public, allocatable, save :: fock_i(:)
    !> Fock basis for the intermediate Hilbert space, valence sector only (length ndim_n_nocore).
    integer(dp), public, allocatable, save :: fock_n(:)
    !> Fock basis for the final Hilbert space (length ndim_f).
    integer(dp), public, allocatable, save :: fock_f(:)

    !> Hamiltonian stored as an array of CSR blocks, one block per MPI column partition.
    type(T_csr), public, allocatable, save :: ham_csr(:)
    !> Transition operator stored as an array of CSR blocks.
    type(T_csr), public, allocatable, save :: tran_csr(:)

    !> Diagonal Lanczos coefficients (alpha) saved between xas_driver/rixs_driver calls.
    real(dp), public, allocatable, save :: krylov_alpha(:)
    !> Off-diagonal Lanczos coefficients (beta) saved between xas_driver/rixs_driver calls.
    real(dp), public, allocatable, save :: krylov_beta(:)

    contains

    !> Allocate hopping_i and initialise all entries to (-1, -1, 0).
    subroutine alloc_hopping_i()
        use m_constants, only: czero
        use m_control,   only: nhopp_i
        implicit none
        integer :: i
        allocate(hopping_i(nhopp_i))
        do i = 1, nhopp_i
            hopping_i(i)%ind1 = -1
            hopping_i(i)%ind2 = -1
            hopping_i(i)%val  = czero
        enddo
        return
    end subroutine alloc_hopping_i

    !> Deallocate hopping_i if allocated.
    subroutine dealloc_hopping_i()
        implicit none
        if (allocated(hopping_i)) deallocate(hopping_i)
        return
    end subroutine dealloc_hopping_i

    !> Allocate coulomb_i and initialise all entries to (-1,-1,-1,-1, 0).
    subroutine alloc_coulomb_i()
        use m_constants, only: czero
        use m_control,   only: ncoul_i
        implicit none
        integer :: i
        allocate(coulomb_i(ncoul_i))
        do i = 1, ncoul_i
            coulomb_i(i)%ind1 = -1
            coulomb_i(i)%ind2 = -1
            coulomb_i(i)%ind3 = -1
            coulomb_i(i)%ind4 = -1
            coulomb_i(i)%val  = czero
        enddo
        return
    end subroutine alloc_coulomb_i

    !> Deallocate coulomb_i if allocated.
    subroutine dealloc_coulomb_i()
        implicit none
        if (allocated(coulomb_i)) deallocate(coulomb_i)
        return
    end subroutine dealloc_coulomb_i

    !> Allocate fock_i (length ndim_i) and fill with -1 as a sentinel.
    subroutine alloc_fock_i()
        use m_control, only: ndim_i
        implicit none
        allocate(fock_i(ndim_i))
        fock_i = -1
        return
    end subroutine alloc_fock_i

    !> Deallocate fock_i if allocated.
    subroutine dealloc_fock_i()
        implicit none
        if (allocated(fock_i)) deallocate(fock_i)
        return
    end subroutine dealloc_fock_i

    !> Allocate hopping_n and initialise all entries to (-1, -1, 0).
    subroutine alloc_hopping_n()
        use m_constants, only: czero
        use m_control,   only: nhopp_n
        implicit none
        integer :: i
        allocate(hopping_n(nhopp_n))
        do i = 1, nhopp_n
            hopping_n(i)%ind1 = -1
            hopping_n(i)%ind2 = -1
            hopping_n(i)%val  = czero
        enddo
        return
    end subroutine alloc_hopping_n

    !> Deallocate hopping_n if allocated.
    subroutine dealloc_hopping_n()
        implicit none
        if (allocated(hopping_n)) deallocate(hopping_n)
        return
    end subroutine dealloc_hopping_n

    !> Allocate coulomb_n and initialise all entries to (-1,-1,-1,-1, 0).
    subroutine alloc_coulomb_n()
        use m_constants, only: czero
        use m_control,   only: ncoul_n
        implicit none
        integer :: i
        allocate(coulomb_n(ncoul_n))
        do i = 1, ncoul_n
            coulomb_n(i)%ind1 = -1
            coulomb_n(i)%ind2 = -1
            coulomb_n(i)%ind3 = -1
            coulomb_n(i)%ind4 = -1
            coulomb_n(i)%val  = czero
        enddo
        return
    end subroutine alloc_coulomb_n

    !> Deallocate coulomb_n if allocated.
    subroutine dealloc_coulomb_n()
        implicit none
        if (allocated(coulomb_n)) deallocate(coulomb_n)
        return
    end subroutine dealloc_coulomb_n

    !> Allocate fock_n (length ndim_n_nocore) and fill with -1.
    subroutine alloc_fock_n()
        use m_control, only: ndim_n_nocore
        implicit none
        allocate(fock_n(ndim_n_nocore))
        fock_n = -1
        return
    end subroutine alloc_fock_n

    !> Deallocate fock_n if allocated.
    subroutine dealloc_fock_n()
        implicit none
        if (allocated(fock_n)) deallocate(fock_n)
        return
    end subroutine dealloc_fock_n

    !> Allocate fock_f (length ndim_f) and fill with -1.
    subroutine alloc_fock_f()
        use m_control, only: ndim_f
        implicit none
        allocate(fock_f(ndim_f))
        fock_f = -1
        return
    end subroutine alloc_fock_f

    !> Deallocate fock_f if allocated.
    subroutine dealloc_fock_f()
        implicit none
        if (allocated(fock_f)) deallocate(fock_f)
        return
    end subroutine dealloc_fock_f

    !> Allocate transop_xas and initialise all entries to (-1, -1, 0).
    subroutine alloc_transop_xas()
        use m_constants, only: czero
        use m_control,   only: ntran_xas
        implicit none
        integer :: i
        allocate(transop_xas(ntran_xas))
        do i = 1, ntran_xas
            transop_xas(i)%ind1 = -1
            transop_xas(i)%ind2 = -1
            transop_xas(i)%val  = czero
        enddo
        return
    end subroutine alloc_transop_xas

    !> Deallocate transop_xas if allocated.
    subroutine dealloc_transop_xas()
        implicit none
        if (allocated(transop_xas)) deallocate(transop_xas)
        return
    end subroutine dealloc_transop_xas

    !> Allocate transop_rixs_i and initialise all entries to (-1, -1, 0).
    subroutine alloc_transop_rixs_i()
        use m_constants, only: czero
        use m_control,   only: ntran_rixs_i
        implicit none
        integer :: i
        allocate(transop_rixs_i(ntran_rixs_i))
        do i = 1, ntran_rixs_i
            transop_rixs_i(i)%ind1 = -1
            transop_rixs_i(i)%ind2 = -1
            transop_rixs_i(i)%val  = czero
        enddo
        return
    end subroutine alloc_transop_rixs_i

    !> Deallocate transop_rixs_i if allocated.
    subroutine dealloc_transop_rixs_i()
        implicit none
        if (allocated(transop_rixs_i)) deallocate(transop_rixs_i)
        return
    end subroutine dealloc_transop_rixs_i

    !> Allocate transop_rixs_f and initialise all entries to (-1, -1, 0).
    subroutine alloc_transop_rixs_f()
        use m_constants, only: czero
        use m_control,   only: ntran_rixs_f
        implicit none
        integer :: i
        allocate(transop_rixs_f(ntran_rixs_f))
        do i = 1, ntran_rixs_f
            transop_rixs_f(i)%ind1 = -1
            transop_rixs_f(i)%ind2 = -1
            transop_rixs_f(i)%val  = czero
        enddo
        return
    end subroutine alloc_transop_rixs_f

    !> Deallocate transop_rixs_f if allocated.
    subroutine dealloc_transop_rixs_f()
        implicit none
        if (allocated(transop_rixs_f)) deallocate(transop_rixs_f)
        return
    end subroutine dealloc_transop_rixs_f

    !> Allocate ham_csr as an array of nblock CSR blocks, each initialised to empty.
    !!
    !! @param[in] nblock  Number of CSR blocks (equals nprocs in the parallel case)
    subroutine alloc_ham_csr(nblock)
        implicit none
        integer, intent(in) :: nblock
        integer :: i
        allocate(ham_csr(nblock))
        do i = 1, nblock
            call init_csr(ham_csr(i))
        enddo
        return
    end subroutine alloc_ham_csr

    !> Deallocate all CSR arrays inside ham_csr and then ham_csr itself.
    !!
    !! @param[in] nblock  Number of CSR blocks
    subroutine dealloc_ham_csr(nblock)
        implicit none
        integer, intent(in) :: nblock
        integer :: i
        do i = 1, nblock
            call dealloc_csr(ham_csr(i))
        enddo
        if (allocated(ham_csr)) deallocate(ham_csr)
        return
    end subroutine dealloc_ham_csr

    !> Allocate tran_csr as an array of nblock CSR blocks, each initialised to empty.
    !!
    !! @param[in] nblock  Number of CSR blocks
    subroutine alloc_tran_csr(nblock)
        implicit none
        integer, intent(in) :: nblock
        integer :: i
        allocate(tran_csr(nblock))
        do i = 1, nblock
            call init_csr(tran_csr(i))
        enddo
        return
    end subroutine alloc_tran_csr

    !> Deallocate all CSR arrays inside tran_csr and then tran_csr itself.
    !!
    !! @param[in] nblock  Number of CSR blocks
    subroutine dealloc_tran_csr(nblock)
        implicit none
        integer, intent(in) :: nblock
        integer :: i
        do i = 1, nblock
            call dealloc_csr(tran_csr(i))
        enddo
        if (allocated(tran_csr)) deallocate(tran_csr)
        return
    end subroutine dealloc_tran_csr

end module m_global
