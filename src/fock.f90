!> Compute the binomial coefficient C(n, m).
!!
!! Uses the multiplicative formula, accumulating the numerator and denominator
!! separately in floating point before rounding to an integer, which avoids
!! overflow for the moderate values of n encountered in ED calculations.
!!
!! @param[in]  n    Total number of spin-orbitals
!! @param[in]  m    Number of electrons
!! @param[out] cnm  C(n, m) = n! / (m! (n-m)!)
subroutine cal_combination(n, m, cnm)
    use m_constants, only: dp

    implicit none

    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(out) :: cnm

    integer :: small
    integer :: large
    integer :: i

    real(dp) :: x
    real(dp) :: y

    small = min(m, n-m)
    large = max(m, n-m)

    x = 1.0_dp
    y = 1.0_dp
    do i = large+1, n
        x = x * dble(i)
    enddo

    do i = 1, small
        y = y * dble(i)
    enddo

    cnm = nint(x / y)

    return
end subroutine cal_combination

!> Apply a fermionic creation or annihilation operator to a Fock state.
!!
!! Each many-body basis state is encoded as a 64-bit integer whose k-th bit
!! (0-indexed) equals 1 if spin-orbital k+1 is occupied.  This routine flips
!! bit (pos-1) and computes the fermionic sign from the Jordan-Wigner string
!! (the parity of the number of occupied orbitals below position pos).
!!
!! The caller is responsible for checking that the operator is valid (i.e.
!! that the target orbital is occupied before annihilation, and empty before
!! creation); no such check is performed here.
!!
!! @param[in]  c_type  Operator type: '+' for creation, '-' for annihilation
!! @param[in]  pos     1-based orbital index on which the operator acts
!! @param[in]  old     Input Fock state (bit-encoded occupation integer)
!! @param[out] new     Output Fock state after applying the operator
!! @param[out] sgn     Fermionic sign: +1 or -1
subroutine make_newfock(c_type, pos, old, new, sgn)
    use m_constants, only: mystd, dp

    implicit none

    character(1), intent(in)  :: c_type
    integer,      intent(in)  :: pos
    integer(dp),  intent(in)  :: old
    integer(dp),  intent(out) :: new
    integer,      intent(out) :: sgn

    integer :: i

    ! Count occupied orbitals below pos to determine the Jordan-Wigner sign
    sgn = 0
    do i = 1, pos-1
        if (btest(old, i-1) .eqv. .true.) sgn = sgn + 1
    enddo
    sgn = mod(sgn, 2)
    sgn = (-1)**sgn

    if (c_type == '+') then
        new = old + (2_dp)**(pos-1)
    elseif (c_type == '-') then
        new = old - (2_dp)**(pos-1)
    else
        write(mystd, "(58a,a)") " fedrixs >>> error of create and destroy operator type: ", c_type
        STOP
    endif

    return
end subroutine make_newfock
