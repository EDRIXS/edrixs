!> Derived types and associated linked-list / CSR sparse-matrix utilities.
!!
!! T_tensor2 and T_tensor4 store the non-zero elements of one-body (hopping)
!! and two-body (Coulomb) operators respectively.  T_col / T_container_col
!! form a singly-linked list used as an intermediate row representation when
!! building the Hamiltonian before converting to T_csr (Compressed Sparse Row).
module m_types
    use m_constants, only: dp, czero
    implicit none

    !> One-body operator element: value at orbital indices (ind1, ind2).
    !! Represents the matrix element t_{ind1,ind2} of a hopping term
    !! c^dagger_{ind2} c_{ind1}.
    type T_tensor2
        integer    :: ind1  !< Row orbital index (annihilation site)
        integer    :: ind2  !< Column orbital index (creation site)
        complex(dp) :: val  !< Matrix element value
    end type T_tensor2

    !> Two-body operator element: value at orbital indices (ind1,ind2,ind3,ind4).
    !! Represents U_{ind1,ind2,ind3,ind4} for the Coulomb term
    !! c^dagger_{ind4} c^dagger_{ind3} c_{ind2} c_{ind1}.
    type T_tensor4
        integer    :: ind1  !< First annihilation orbital index
        integer    :: ind2  !< Second annihilation orbital index
        integer    :: ind3  !< First creation orbital index
        integer    :: ind4  !< Second creation orbital index
        complex(dp) :: val  !< Matrix element value
    end type T_tensor4

    !> A single node in a per-row linked list of sparse matrix entries.
    !! Used as an intermediate structure when constructing rows of the
    !! Hamiltonian before the final conversion to CSR format.
    type T_col
        integer    :: col              !< Column index of this non-zero entry
        complex(dp) :: val             !< Value of this non-zero entry
        type(T_col), pointer :: next   !< Pointer to the next column entry in this row
    end type T_col

    !> Container holding a pointer to the head T_col node for one row.
    !! Needed because Fortran arrays cannot directly hold polymorphic pointers,
    !! so this wrapper type enables an array of per-row linked-list heads.
    type T_container_col
        type(T_col), pointer :: col_ptr  !< Pointer to the first T_col node in the row
    end type T_container_col

    !> Compressed Sparse Row (CSR) representation of a complex sparse matrix.
    !!
    !! The matrix is stored with standard CSR arrays (iaa, jaa, aa) plus
    !! row_shift and col_shift offsets to support distributed storage where
    !! each MPI rank owns a contiguous block of rows and columns.
    type T_csr
        integer :: nnz        !< Number of non-zero elements
        integer :: m          !< Number of rows owned by this rank
        integer :: row_shift  !< Global row index of the first local row minus 1
        integer :: col_shift  !< Global column index of the first local column minus 1
        integer,    allocatable :: iaa(:)  !< Row pointer array (length m+1); iaa(i) is the index in jaa/aa of the first entry of row i
        integer,    allocatable :: jaa(:)  !< Column index array (length nnz)
        complex(dp), allocatable :: aa(:)  !< Non-zero value array (length nnz)
    end type T_csr

    contains

    !> Allocate and initialise a new T_col node.
    !!
    !! @param[out] self  Pointer that will point to the newly allocated node
    !! @param[in]  col   Column index for this node
    !! @param[in]  val   Value to store in this node
    subroutine init_col(self, col, val)
        implicit none
        type(T_col), pointer    :: self
        integer,     intent(in) :: col
        complex(dp), intent(in) :: val

        allocate(self)
        self%col  = col
        self%val  = val
        self%next => null()

        return
    end subroutine init_col

    !> Insert a (col, val) entry into a sorted linked list representing one row.
    !!
    !! If a node with the same column already exists its value is accumulated
    !! (+=).  Otherwise a new node is inserted in ascending column order so
    !! the list remains sorted, which simplifies the subsequent CSR conversion.
    !!
    !! @param[inout] self  Head pointer of the linked list for this row
    !! @param[in]    col   Column index to insert
    !! @param[in]    val   Value to insert or accumulate
    subroutine insert_into_row(self, col, val)
        implicit none
        type(T_col), pointer    :: self
        integer,     intent(in) :: col
        complex(dp), intent(in) :: val

        type(T_col), pointer :: curr
        type(T_col), pointer :: prev
        type(T_col), pointer :: tmp

        curr => self
        prev => null()
        do while (associated(curr))
            if (col == curr%col) then
                curr%val = curr%val + val
                return
            endif
            if (col < curr%col) then
                allocate(tmp)
                tmp%col  = col
                tmp%val  = val
                tmp%next => curr
                if (.not. associated(prev)) then
                    self => tmp
                    return
                endif
                prev%next => tmp
                return
            endif
            prev => curr
            curr => curr%next
        enddo

        allocate(tmp)
        tmp%col  = col
        tmp%val  = val
        tmp%next => null()
        prev%next => tmp

        return
    end subroutine insert_into_row

    !> Initialise a T_csr struct to empty (zero sizes, zero shifts).
    !!
    !! @param[inout] self  The T_csr object to initialise
    subroutine init_csr(self)
        implicit none
        type(T_csr) :: self

        self%m         = 0
        self%nnz       = 0
        self%row_shift = 0
        self%col_shift = 0

        return
    end subroutine init_csr

    !> Allocate the CSR arrays (iaa, jaa, aa) for a matrix of given size.
    !!
    !! @param[in]    m     Number of rows
    !! @param[in]    nnz   Number of non-zero elements
    !! @param[inout] self  T_csr object whose arrays are allocated
    subroutine alloc_csr(m, nnz, self)
        implicit none
        integer,    intent(in) :: m
        integer,    intent(in) :: nnz
        type(T_csr)            :: self

        integer :: istat

        allocate(self%iaa(m+1), stat=istat)
        allocate(self%jaa(nnz), stat=istat)
        allocate(self%aa(nnz),  stat=istat)

        self%iaa = 0
        self%jaa = 0
        self%aa  = czero

        return
    end subroutine alloc_csr

    !> Deallocate the CSR arrays of a T_csr object if they are allocated.
    !!
    !! @param[inout] self  T_csr object to free
    subroutine dealloc_csr(self)
        implicit none
        type(T_csr) :: self

        if (allocated(self%iaa)) deallocate(self%iaa)
        if (allocated(self%jaa)) deallocate(self%jaa)
        if (allocated(self%aa))  deallocate(self%aa)

        return
    end subroutine dealloc_csr

    !> Copy non-zero entries from a T_csr block into a dense matrix.
    !!
    !! Only the entries owned by m_csr are written; the caller must ensure
    !! m_full is large enough and pre-initialised if a full gather is needed.
    !!
    !! @param[in]  m       Number of rows of the full matrix
    !! @param[in]  n       Number of columns of the full matrix
    !! @param[in]  m_csr   Source CSR block (may be a sub-block of the full matrix)
    !! @param[out] m_full  Destination dense matrix (m x n)
    subroutine csr_to_full(m, n, m_csr, m_full)
        implicit none
        integer,     intent(in)  :: m
        integer,     intent(in)  :: n
        type(T_csr)              :: m_csr
        complex(dp), intent(out) :: m_full(m, n)

        integer :: i, j, begin_col, end_col

        do i = 1, m_csr%m
            begin_col = m_csr%iaa(i)
            end_col   = m_csr%iaa(i+1) - 1
            if (begin_col > end_col) cycle
            do j = begin_col, end_col
                m_full(i + m_csr%row_shift, m_csr%jaa(j)) = m_csr%aa(j)
            enddo
        enddo

        return
    end subroutine csr_to_full

end module m_types
