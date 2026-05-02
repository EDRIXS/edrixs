!> Build the many-body Hamiltonian for the initial (no core-hole) configuration.
!!
!! Iterates over each Fock basis state owned by this MPI rank and applies every
!! non-zero hopping and Coulomb operator.  For each operator the resulting
!! state is located in the Fock basis with binary_search; if found, the matrix
!! element (with Jordan-Wigner sign from make_newfock) is inserted into a
!! per-row linked list.  The loop runs twice: the first pass counts non-zeros
!! per row so CSR storage can be allocated exactly; the second pass fills the
!! values.
!!
!! The optional diagonal shift omega (added only when ham_sign < 0) supports
!! building the resolvent (H - omega I) needed for RIXS intermediate states.
!! ham_sign = +1 builds H; ham_sign = -1 builds -(H - omega I).
!!
!! Parallelism: rows are distributed across nprocs ranks according to
!! end_indx(1,1,myid+1)..end_indx(2,1,myid+1).  Column blocks are distributed
!! across nblock blocks aligned with the same partition.
!!
!! @param[in]  ncfgs     Total dimension of the Hilbert space
!! @param[in]  fock      Fock basis array (length ncfgs), sorted ascending
!! @param[in]  nblock    Number of CSR column blocks (equals nprocs in the parallel case)
!! @param[in]  end_indx  Global index ranges per rank: end_indx(start/end, row/col, rank)
!! @param[in]  nhopp     Number of non-zero hopping terms
!! @param[in]  hopping   Hopping operator elements (length nhopp)
!! @param[in]  ncoul     Number of non-zero Coulomb terms
!! @param[in]  coulomb   Coulomb operator elements (length ncoul)
!! @param[in]  omega     Diagonal energy shift (used when ham_sign < 0)
!! @param[in]  ham_sign  +1.0 to build H; -1.0 to build -(H - omega I)
!! @param[out] ham_csr   Output Hamiltonian as an array of nblock CSR blocks
subroutine build_ham_i(ncfgs, fock, nblock, end_indx, nhopp, hopping, &
                       ncoul, coulomb, omega, ham_sign, ham_csr)
    use m_constants, only: dp
    use m_control,   only: myid, master, nprocs
    use m_types
    use mpi

    implicit none

    integer,          intent(in) :: ncfgs
    integer(dp),      intent(in) :: fock(ncfgs)
    integer,          intent(in) :: nblock
    integer,          intent(in) :: end_indx(2, 2, nprocs)
    integer,          intent(in) :: nhopp
    type(T_tensor2),  intent(in) :: hopping(nhopp)
    integer,          intent(in) :: ncoul
    type(T_tensor4),  intent(in) :: coulomb(ncoul)
    complex(dp),      intent(in) :: omega
    real(dp),         intent(in) :: ham_sign
    type(T_csr)                  :: ham_csr(nblock)

    integer, external :: binary_search
    integer(dp) :: old, new
    integer :: icfg, jcfg
    integer :: alpha, beta, gama, delta
    integer :: tot_sign, sgn
    integer :: i, j
    integer :: which_block, row_shift, begin_col

    integer, allocatable :: num_of_cols(:,:)
    integer :: num_of_tot(nblock)

    type(T_col), pointer :: next_col, curr_col
    type(T_container_col) :: cols(nblock)

    which_block = 0
    row_shift   = end_indx(1,1,myid+1) - 1
    allocate(num_of_cols(end_indx(2,1,myid+1)-end_indx(1,1,myid+1)+1, nblock))
    num_of_cols = 0
    num_of_tot  = 0

    ! --- First pass: count non-zeros per row ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo

        do i = 1, nhopp
            alpha = hopping(i)%ind1
            beta  = hopping(i)%ind2
            if (.not. btest(fock(icfg), alpha-1)) cycle
            old = fock(icfg)
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn
            old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            jcfg = binary_search(ncfgs, fock, new)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            endif
        enddo

        do i = 1, ncoul
            alpha = coulomb(i)%ind1
            beta  = coulomb(i)%ind2
            gama  = coulomb(i)%ind3
            delta = coulomb(i)%ind4
            if (.not. btest(fock(icfg), alpha-1)) cycle
            old = fock(icfg)
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (.not. btest(old, beta-1)) cycle
            call make_newfock('-', beta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, gama-1)) cycle
            call make_newfock('+', gama, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, delta-1)) cycle
            call make_newfock('+', delta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            jcfg = binary_search(ncfgs, fock, new)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            endif
        enddo

        ! Diagonal shift omega (for resolvent H - omega I)
        if (ham_sign < 0) then
            jcfg = icfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, omega)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, omega)
            endif
        endif

        do i = 1, nblock
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                num_of_cols(icfg-row_shift, i) = num_of_cols(icfg-row_shift, i) + 1
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
            num_of_tot(i) = num_of_tot(i) + num_of_cols(icfg-row_shift, i)
        enddo
    enddo

    if (myid == master) print *, "    Allocate memory for ham_csr ..."
    do i = 1, nblock
        if (num_of_tot(i) /= 0) then
            ham_csr(i)%m         = end_indx(2,1,myid+1) - end_indx(1,1,myid+1) + 1
            ham_csr(i)%nnz       = num_of_tot(i)
            ham_csr(i)%row_shift = end_indx(1,1,myid+1) - 1
            ham_csr(i)%col_shift = end_indx(1,2,i) - 1
            call alloc_csr(ham_csr(i)%m, ham_csr(i)%nnz, ham_csr(i))
            ham_csr(i)%iaa(1) = 1
            do j = 2, ham_csr(i)%m+1
                ham_csr(i)%iaa(j) = ham_csr(i)%iaa(j-1) + num_of_cols(j-1, i)
            enddo
        endif
    enddo
    deallocate(num_of_cols)

    ! --- Second pass: fill CSR values ---
    if (myid == master) print *, "    Really building Hamiltonian begin here ..."
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo

        do i = 1, nhopp
            alpha = hopping(i)%ind1
            beta  = hopping(i)%ind2
            if (.not. btest(fock(icfg), alpha-1)) cycle
            old = fock(icfg)
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            jcfg = binary_search(ncfgs, fock, new)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            endif
        enddo

        do i = 1, ncoul
            alpha = coulomb(i)%ind1
            beta  = coulomb(i)%ind2
            gama  = coulomb(i)%ind3
            delta = coulomb(i)%ind4
            if (.not. btest(fock(icfg), alpha-1)) cycle
            old = fock(icfg)
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (.not. btest(old, beta-1)) cycle
            call make_newfock('-', beta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, gama-1)) cycle
            call make_newfock('+', gama, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, delta-1)) cycle
            call make_newfock('+', delta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            jcfg = binary_search(ncfgs, fock, new)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            endif
        enddo

        if (ham_sign < 0) then
            jcfg = icfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, omega)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, omega)
            endif
        endif

        do i = 1, nblock
            if (num_of_tot(i) == 0) cycle
            begin_col = ham_csr(i)%iaa(icfg-row_shift) - 1
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                begin_col = begin_col + 1
                ham_csr(i)%jaa(begin_col) = next_col%col
                ham_csr(i)%aa(begin_col)  = next_col%val
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
        enddo
    enddo

    return
end subroutine build_ham_i

!> Build the many-body Hamiltonian for the intermediate (core-hole) configuration.
!!
!! Extends build_ham_i to the product Hilbert space
!! H_n = H_val-only (x) H_core-hole, where H_core-hole has exactly one core
!! orbital empty.  The num_core_orbs possible core-hole configurations are
!! enumerated via fock_core(1..num_core_orbs), each encoding one missing
!! core orbital as a bit pattern.  The composite Fock index is linearised as
!!   icfg = (core_indx - 1) * ncfgs + val_indx
!! so the same CSR row-block structure as build_ham_i applies.
!!
!! @param[in]  ncfgs         Total dimension of the valence-sector Hilbert space
!! @param[in]  fock          Valence-sector Fock basis (length ncfgs)
!! @param[in]  num_val_orbs  Number of valence spin-orbitals
!! @param[in]  num_core_orbs Number of core spin-orbitals
!! @param[in]  nblock        Number of CSR column blocks
!! @param[in]  end_indx      Global index ranges per rank
!! @param[in]  nhopp         Number of non-zero hopping terms
!! @param[in]  hopping       Hopping operator elements
!! @param[in]  ncoul         Number of non-zero Coulomb terms
!! @param[in]  coulomb       Coulomb operator elements
!! @param[in]  omega         Diagonal energy shift (used when ham_sign < 0)
!! @param[in]  ham_sign      +1.0 to build H; -1.0 to build -(H - omega I)
!! @param[out] ham_csr       Output Hamiltonian as an array of nblock CSR blocks
subroutine build_ham_n(ncfgs, fock, num_val_orbs, num_core_orbs, nblock, &
       end_indx, nhopp, hopping, ncoul, coulomb, omega, ham_sign, ham_csr)
    use m_constants, only: dp
    use m_control,   only: myid, master, nprocs
    use m_types
    use mpi

    implicit none

    integer,         intent(in) :: ncfgs
    integer(dp),     intent(in) :: fock(ncfgs)
    integer,         intent(in) :: num_val_orbs
    integer,         intent(in) :: num_core_orbs
    integer,         intent(in) :: nblock
    integer,         intent(in) :: end_indx(2, 2, nprocs)
    integer,         intent(in) :: nhopp
    type(T_tensor2), intent(in) :: hopping(nhopp)
    integer,         intent(in) :: ncoul
    type(T_tensor4), intent(in) :: coulomb(ncoul)
    complex(dp),     intent(in) :: omega
    real(dp),        intent(in) :: ham_sign
    type(T_csr)                 :: ham_csr(nblock)

    integer, external :: binary_search
    integer(dp) :: old, new
    integer :: icfg, jcfg, kcfg, lcfg
    integer :: alpha, beta, gama, delta
    integer :: tot_sign, sgn
    integer :: i, j
    integer :: which_block, row_shift, begin_col
    integer :: core_indx, val_indx
    integer(dp) :: fock_core(num_core_orbs)
    integer(dp) :: new_core, new_val

    integer, allocatable :: num_of_cols(:,:)
    integer :: num_of_tot(nblock)

    type(T_col), pointer :: next_col, curr_col
    type(T_container_col) :: cols(nblock)

    which_block = 0
    row_shift   = end_indx(1,1,myid+1) - 1
    allocate(num_of_cols(end_indx(2,1,myid+1)-end_indx(1,1,myid+1)+1, nblock))
    num_of_cols = 0
    num_of_tot  = 0

    ! Build the num_core_orbs possible single-core-hole Fock states.
    ! fock_core(k) has all core bits set except bit (num_core_orbs - k).
    do i = 1, num_core_orbs
        fock_core(i) = (2_dp)**num_core_orbs - 1_dp - (2_dp)**(num_core_orbs-i)
    enddo

    ! --- First pass: count non-zeros per row ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        core_indx = (icfg-1)/ncfgs + 1
        val_indx  = mod(icfg-1, ncfgs) + 1
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo

        do i = 1, nhopp
            alpha = hopping(i)%ind1
            beta  = hopping(i)%ind2
            ! Composite state: core bits shifted above valence bits
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            endif
        enddo

        do i = 1, ncoul
            alpha = coulomb(i)%ind1
            beta  = coulomb(i)%ind2
            gama  = coulomb(i)%ind3
            delta = coulomb(i)%ind4
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (.not. btest(old, beta-1)) cycle
            call make_newfock('-', beta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, gama-1)) cycle
            call make_newfock('+', gama, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, delta-1)) cycle
            call make_newfock('+', delta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            endif
        enddo

        if (ham_sign < 0) then
            jcfg = icfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, omega)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, omega)
            endif
        endif

        do i = 1, nblock
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                num_of_cols(icfg-row_shift, i) = num_of_cols(icfg-row_shift, i) + 1
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
            num_of_tot(i) = num_of_tot(i) + num_of_cols(icfg-row_shift, i)
        enddo
    enddo

    if (myid == master) print *, "    Allocate memory for ham_csr ..."
    do i = 1, nblock
        if (num_of_tot(i) /= 0) then
            ham_csr(i)%m         = end_indx(2,1,myid+1) - end_indx(1,1,myid+1) + 1
            ham_csr(i)%nnz       = num_of_tot(i)
            ham_csr(i)%row_shift = end_indx(1,1,myid+1) - 1
            ham_csr(i)%col_shift = end_indx(1,2,i) - 1
            call alloc_csr(ham_csr(i)%m, ham_csr(i)%nnz, ham_csr(i))
            ham_csr(i)%iaa(1) = 1
            do j = 2, ham_csr(i)%m+1
                ham_csr(i)%iaa(j) = ham_csr(i)%iaa(j-1) + num_of_cols(j-1, i)
            enddo
        endif
    enddo
    deallocate(num_of_cols)

    ! --- Second pass: fill CSR values ---
    if (myid == master) print *, "    Really building Hamiltonian begin here ..."
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        core_indx = (icfg-1)/ncfgs + 1
        val_indx  = mod(icfg-1, ncfgs) + 1
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo

        do i = 1, nhopp
            alpha = hopping(i)%ind1
            beta  = hopping(i)%ind2
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, hopping(i)%val * tot_sign * ham_sign)
            endif
        enddo

        do i = 1, ncoul
            alpha = coulomb(i)%ind1
            beta  = coulomb(i)%ind2
            gama  = coulomb(i)%ind3
            delta = coulomb(i)%ind4
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (.not. btest(old, beta-1)) cycle
            call make_newfock('-', beta, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, gama-1)) cycle
            call make_newfock('+', gama, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, delta-1)) cycle
            call make_newfock('+', delta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, coulomb(i)%val * tot_sign * ham_sign)
            endif
        enddo

        if (ham_sign < 0) then
            jcfg = icfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, omega)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, omega)
            endif
        endif

        do i = 1, nblock
            if (num_of_tot(i) == 0) cycle
            begin_col = ham_csr(i)%iaa(icfg-row_shift) - 1
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                begin_col = begin_col + 1
                ham_csr(i)%jaa(begin_col) = next_col%col
                ham_csr(i)%aa(begin_col)  = next_col%val
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
        enddo
    enddo

    return
end subroutine build_ham_n

!> Build the transition operator from the initial space to the intermediate space.
!!
!! Constructs the sparse matrix representation of the photon absorption operator
!! T: H_i -> H_n, where H_i is the initial (no core-hole) space and H_n is
!! the intermediate (core-hole) space.  The composite H_n index interleaves
!! core-hole configurations with the valence-sector index.
!!
!! The fock_left array contains valence-sector states from H_n (rows) and
!! fock_right contains the initial states from H_i (columns).  The routine
!! uses the same two-pass strategy as build_ham_i.
!!
!! @param[in]  mcfgs         Dimension of the left (intermediate valence) Hilbert space
!! @param[in]  ncfgs         Dimension of the right (initial) Hilbert space
!! @param[in]  fock_left     Fock basis of the intermediate valence sector (length mcfgs)
!! @param[in]  fock_right    Fock basis of the initial sector (length ncfgs)
!! @param[in]  num_val_orbs  Number of valence spin-orbitals
!! @param[in]  num_core_orbs Number of core spin-orbitals
!! @param[in]  nblock        Number of CSR column blocks
!! @param[in]  end_indx      Global index ranges per rank
!! @param[in]  ntran         Number of non-zero transition operator elements
!! @param[in]  transop       Transition operator elements (length ntran)
!! @param[out] tran_csr      Output transition operator as an array of nblock CSR blocks
subroutine build_transop_i(mcfgs, ncfgs, fock_left, fock_right, num_val_orbs, &
                     num_core_orbs, nblock, end_indx, ntran, transop, tran_csr)
    use m_constants, only: dp
    use mpi
    use m_control,   only: myid, nprocs
    use m_types

    implicit none

    integer,         intent(in) :: mcfgs
    integer,         intent(in) :: ncfgs
    integer(dp),     intent(in) :: fock_left(mcfgs)
    integer(dp),     intent(in) :: fock_right(ncfgs)
    integer,         intent(in) :: num_val_orbs
    integer,         intent(in) :: num_core_orbs
    integer,         intent(in) :: nblock
    integer,         intent(in) :: end_indx(2, 2, nprocs)
    integer,         intent(in) :: ntran
    type(T_tensor2), intent(in) :: transop(ntran)
    type(T_csr)                 :: tran_csr(nblock)

    integer, external :: binary_search
    integer :: icfg, jcfg
    integer :: alpha, beta
    integer(dp) :: old, new
    integer :: tot_sign, sgn
    integer :: i, j
    integer :: which_block, row_shift, begin_col
    integer :: core_indx, val_indx
    integer(dp) :: fock_core(num_core_orbs)
    integer(dp) :: new_core, new_val

    integer, allocatable :: num_of_cols(:,:)
    integer :: num_of_tot(nblock)

    type(T_col), pointer :: next_col, curr_col
    type(T_container_col) :: cols(nblock)

    which_block = 0
    row_shift   = end_indx(1,1,myid+1) - 1
    allocate(num_of_cols(end_indx(2,1,myid+1)-end_indx(1,1,myid+1)+1, nblock))
    num_of_cols = 0
    num_of_tot  = 0

    do i = 1, num_core_orbs
        fock_core(i) = (2_dp)**num_core_orbs - (1_dp) - (2_dp)**(num_core_orbs-i)
    enddo

    ! --- First pass: count non-zeros ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        core_indx = (icfg-1)/mcfgs + 1
        val_indx  = mod(icfg-1, mcfgs) + 1
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo
        do i = 1, ntran
            alpha = transop(i)%ind1
            beta  = transop(i)%ind2
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock_left(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            jcfg = binary_search(ncfgs, fock_right, new_val)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            endif
        enddo

        do i = 1, nblock
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                num_of_cols(icfg-row_shift, i) = num_of_cols(icfg-row_shift, i) + 1
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
            num_of_tot(i) = num_of_tot(i) + num_of_cols(icfg-row_shift, i)
        enddo
    enddo

    do i = 1, nblock
        if (num_of_tot(i) /= 0) then
            tran_csr(i)%m         = end_indx(2,1,myid+1) - end_indx(1,1,myid+1) + 1
            tran_csr(i)%nnz       = num_of_tot(i)
            tran_csr(i)%row_shift = end_indx(1,1,myid+1) - 1
            tran_csr(i)%col_shift = end_indx(1,2,i) - 1
            call alloc_csr(tran_csr(i)%m, tran_csr(i)%nnz, tran_csr(i))
            tran_csr(i)%iaa(1) = 1
            do j = 2, tran_csr(i)%m+1
                tran_csr(i)%iaa(j) = tran_csr(i)%iaa(j-1) + num_of_cols(j-1, i)
            enddo
        endif
    enddo
    deallocate(num_of_cols)

    ! --- Second pass: fill values ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        core_indx = (icfg-1)/mcfgs + 1
        val_indx  = mod(icfg-1, mcfgs) + 1
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo
        do i = 1, ntran
            alpha = transop(i)%ind1
            beta  = transop(i)%ind2
            old = fock_core(core_indx) * (2_dp)**num_val_orbs + fock_left(val_indx)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            jcfg = binary_search(ncfgs, fock_right, new_val)
            if (jcfg == -1) cycle
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            endif
        enddo

        do i = 1, nblock
            if (num_of_tot(i) == 0) cycle
            begin_col = tran_csr(i)%iaa(icfg-row_shift) - 1
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                begin_col = begin_col + 1
                tran_csr(i)%jaa(begin_col) = next_col%col
                tran_csr(i)%aa(begin_col)  = next_col%val
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
        enddo
    enddo

    return
end subroutine build_transop_i

!> Build the transition operator from the intermediate space to the final space.
!!
!! Constructs the photon emission operator T: H_n -> H_f.  The intermediate
!! space H_n has all core orbitals fully occupied (full core), so the
!! composite Fock state for the row is built with ((2^num_core_orbs - 1) << num_val_orbs)
!! ORed with the valence part from fock_left.  The column states in H_f are
!! found by searching fock_right after decomposing the result into core and
!! valence sectors.
!!
!! @param[in]  mcfgs         Dimension of the left (final) Hilbert space
!! @param[in]  ncfgs         Dimension of the right (intermediate valence) Hilbert space
!! @param[in]  fock_left     Fock basis of the final sector (length mcfgs)
!! @param[in]  fock_right    Fock basis of the intermediate valence sector (length ncfgs)
!! @param[in]  num_val_orbs  Number of valence spin-orbitals
!! @param[in]  num_core_orbs Number of core spin-orbitals
!! @param[in]  nblock        Number of CSR column blocks
!! @param[in]  end_indx      Global index ranges per rank
!! @param[in]  ntran         Number of non-zero transition operator elements
!! @param[in]  transop       Transition operator elements (length ntran)
!! @param[out] tran_csr      Output transition operator as an array of nblock CSR blocks
subroutine build_transop_f(mcfgs, ncfgs, fock_left, fock_right, num_val_orbs, &
                     num_core_orbs, nblock, end_indx, ntran, transop, tran_csr)
    use m_constants, only: dp
    use mpi
    use m_control,   only: myid, nprocs
    use m_types

    implicit none

    integer,         intent(in) :: mcfgs
    integer,         intent(in) :: ncfgs
    integer(dp),     intent(in) :: fock_left(mcfgs)
    integer(dp),     intent(in) :: fock_right(ncfgs)
    integer,         intent(in) :: num_val_orbs
    integer,         intent(in) :: num_core_orbs
    integer,         intent(in) :: nblock
    integer,         intent(in) :: end_indx(2, 2, nprocs)
    integer,         intent(in) :: ntran
    type(T_tensor2), intent(in) :: transop(ntran)
    type(T_csr)                 :: tran_csr(nblock)

    integer, external :: binary_search
    integer :: icfg, jcfg, kcfg, lcfg
    integer :: alpha, beta
    integer(dp) :: old, new
    integer :: tot_sign, sgn
    integer :: i, j
    integer :: which_block, row_shift, begin_col
    integer :: core_indx, val_indx
    integer(dp) :: fock_core(num_core_orbs)
    integer(dp) :: new_core, new_val

    integer, allocatable :: num_of_cols(:,:)
    integer :: num_of_tot(nblock)

    type(T_col), pointer :: next_col, curr_col
    type(T_container_col) :: cols(nblock)

    which_block = 0
    row_shift   = end_indx(1,1,myid+1) - 1
    allocate(num_of_cols(end_indx(2,1,myid+1)-end_indx(1,1,myid+1)+1, nblock))
    num_of_cols = 0
    num_of_tot  = 0

    do i = 1, num_core_orbs
        fock_core(i) = (2_dp)**num_core_orbs - (1_dp) - (2_dp)**(num_core_orbs-i)
    enddo

    ! --- First pass: count non-zeros ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo
        do i = 1, ntran
            alpha = transop(i)%ind1
            beta  = transop(i)%ind2
            ! Full core: all core bits set
            old = ((2_dp)**num_core_orbs - (1_dp)) * (2_dp)**num_val_orbs + fock_left(icfg)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock_right, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            endif
        enddo

        do i = 1, nblock
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                num_of_cols(icfg-row_shift, i) = num_of_cols(icfg-row_shift, i) + 1
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
            num_of_tot(i) = num_of_tot(i) + num_of_cols(icfg-row_shift, i)
        enddo
    enddo

    do i = 1, nblock
        if (num_of_tot(i) /= 0) then
            tran_csr(i)%m         = end_indx(2,1,myid+1) - end_indx(1,1,myid+1) + 1
            tran_csr(i)%nnz       = num_of_tot(i)
            tran_csr(i)%row_shift = end_indx(1,1,myid+1) - 1
            tran_csr(i)%col_shift = end_indx(1,2,i) - 1
            call alloc_csr(tran_csr(i)%m, tran_csr(i)%nnz, tran_csr(i))
            tran_csr(i)%iaa(1) = 1
            do j = 2, tran_csr(i)%m+1
                tran_csr(i)%iaa(j) = tran_csr(i)%iaa(j-1) + num_of_cols(j-1, i)
            enddo
        endif
    enddo
    deallocate(num_of_cols)

    ! --- Second pass: fill values ---
    do icfg = end_indx(1,1,myid+1), end_indx(2,1,myid+1)
        core_indx = (icfg-1)/mcfgs + 1
        val_indx  = mod(icfg-1, mcfgs) + 1
        do i = 1, nblock
            cols(i)%col_ptr => null()
        enddo
        do i = 1, ntran
            alpha = transop(i)%ind1
            beta  = transop(i)%ind2
            old = ((2_dp)**num_core_orbs - (1_dp)) * (2_dp)**num_val_orbs + fock_left(icfg)
            if (.not. btest(old, alpha-1)) cycle
            tot_sign = 1
            call make_newfock('-', alpha, old, new, sgn)
            tot_sign = tot_sign * sgn; old = new
            if (btest(old, beta-1)) cycle
            call make_newfock('+', beta, old, new, sgn)
            tot_sign = tot_sign * sgn
            new_core = new / ((2_dp)**num_val_orbs)
            new_val  = new - new_core * (2_dp)**num_val_orbs
            kcfg = binary_search(num_core_orbs, fock_core, new_core)
            lcfg = binary_search(ncfgs, fock_right, new_val)
            if (kcfg == -1 .or. lcfg == -1) cycle
            jcfg = (kcfg-1)*ncfgs + lcfg
            if (nblock == nprocs) then
                do j = 1, nprocs
                    if (jcfg >= end_indx(1,2,j) .and. jcfg <= end_indx(2,2,j)) which_block = j
                enddo
            elseif (nblock == 1) then
                which_block = 1
            else
                print *, "nblock /= 1 and nblock /= nprocs"; STOP
            endif
            if (.not. associated(cols(which_block)%col_ptr)) then
                call init_col(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            else
                call insert_into_row(cols(which_block)%col_ptr, jcfg, transop(i)%val * tot_sign)
            endif
        enddo

        do i = 1, nblock
            if (num_of_tot(i) == 0) cycle
            begin_col = tran_csr(i)%iaa(icfg-row_shift) - 1
            next_col => cols(i)%col_ptr
            do while (associated(next_col))
                begin_col = begin_col + 1
                tran_csr(i)%jaa(begin_col) = next_col%col
                tran_csr(i)%aa(begin_col)  = next_col%val
                curr_col => next_col
                next_col => next_col%next
                deallocate(curr_col)
            enddo
        enddo
    enddo

    return
end subroutine build_transop_f
