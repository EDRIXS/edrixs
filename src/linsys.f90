!> Parallel MINRES solver for the complex symmetric linear system A x = b.
!!
!! Solves (H - omega) x = b where H is the intermediate-state Hamiltonian
!! and omega = omega_in + i*gamma_in encodes the incident photon energy and
!! core-hole lifetime broadening.  The system is complex symmetric (not
!! Hermitian) because omega has a non-zero imaginary part.
!!
!! The algorithm is the Minimum Residual (MINRES) method adapted for
!! complex systems, using pspmv_csr for the distributed matrix-vector
!! product.  Iteration stops when the relative residual drops below
!! linsys_tol or linsys_max iterations are reached.
!!
!! On exit the actual residual |H x - b| is printed by the master rank
!! as a diagnostic check.
!!
!! @param[in]  nblock    Number of CSR column blocks
!! @param[in]  end_indx  Global index ranges per rank: end_indx(start/end, row/col, rank)
!! @param[in]  needed    Communication map from get_needed_indx
!! @param[in]  nloc      Number of local rows / columns
!! @param[in]  ham       Local CSR blocks of the shifted Hamiltonian (H - omega I)
!! @param[in]  b_vec     Right-hand-side vector (length nloc)
!! @param[out] x_vec     Solution vector (length nloc)
!! @param[out] info      Return code: 0 = converged, 1 = max iterations reached
subroutine pminres_csr(nblock, end_indx, needed, nloc, ham, b_vec, x_vec, info)
    use m_constants, only: dp, zero, czero, cone, mystd
    use m_control,   only: linsys_tol, linsys_max, nprocs, myid, master, new_comm
    use m_types
    use mpi

    implicit none

    integer, intent(in)       :: nblock
    integer, intent(in)       :: end_indx(2, 2, nprocs)
    integer, intent(in)       :: needed(nprocs, nprocs)
    integer, intent(in)       :: nloc
    type(T_csr)               :: ham(nblock)
    complex(dp), intent(in)   :: b_vec(nloc)
    complex(dp), intent(out)  :: x_vec(nloc)
    integer,     intent(out)  :: info

    integer  :: curr, prev, tmp_indx, ierror, iter
    real(dp) :: beta0, beta1, beta0_mpi, beta1_mpi
    complex(dp) :: sigma0, sigma1, gamma0, gamma1
    complex(dp) :: eta, delta, res
    complex(dp) :: alpha, alpha_mpi
    complex(dp) :: rho1, rho2, rho3

    complex(dp), allocatable :: v_vec(:,:)
    complex(dp), allocatable :: w_vec(:,:)
    complex(dp), allocatable :: tmp_vec(:)

    info = 0

    allocate(v_vec(nloc, 2))
    allocate(w_vec(nloc, 2))
    allocate(tmp_vec(nloc))

    x_vec   = czero
    v_vec   = czero
    w_vec   = czero
    tmp_vec = czero
    sigma0  = czero
    sigma1  = czero
    gamma0  = cone
    gamma1  = cone
    beta0   = zero
    beta1   = zero
    beta0_mpi = zero
    beta1_mpi = zero
    alpha   = czero
    alpha_mpi = czero
    rho1    = czero
    rho2    = czero
    rho3    = czero
    delta   = czero

    curr = 1
    prev = 2

    v_vec(:, curr) = b_vec
    beta0 = real(dot_product(v_vec(:, curr), v_vec(:, curr)))
    call MPI_ALLREDUCE(beta0, beta0_mpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, new_comm, ierror)
    beta0 = sqrt(beta0_mpi)
    eta   = beta0
    res   = beta0

    iter = 0
    do while (abs(res) > linsys_tol .and. iter < linsys_max)
        iter = iter + 1
        if (myid == master .and. mod(iter, 50) == 0) then
            write(mystd, "(a25, i5, E10.2, a5, E10.2)") "PMINRES iteration:  ", iter, abs(res), "-->", linsys_tol
        endif
        v_vec(:, curr) = v_vec(:, curr) / beta0
        tmp_vec = czero
        call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, v_vec(:, curr), tmp_vec)
        alpha = dot_product(v_vec(:, curr), tmp_vec)
        call MPI_ALLREDUCE(alpha, alpha_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
        alpha = alpha_mpi
        v_vec(:, prev) = tmp_vec - alpha * v_vec(:, curr) - beta0 * v_vec(:, prev)
        beta1 = real(dot_product(v_vec(:, prev), v_vec(:, prev)))
        call MPI_ALLREDUCE(beta1, beta1_mpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, new_comm, ierror)
        beta1 = sqrt(beta1_mpi)

        delta = gamma0 * alpha - gamma1 * sigma0 * beta0
        rho1  = sqrt(delta*delta + beta1*beta1)
        rho2  = sigma0 * alpha + gamma1 * gamma0 * beta0
        rho3  = sigma1 * beta0

        gamma1 = gamma0
        sigma1 = sigma0
        gamma0 = delta / rho1
        sigma0 = beta1 / rho1

        w_vec(:, prev) = (v_vec(:, curr) - rho3 * w_vec(:, prev) - rho2 * w_vec(:, curr)) / rho1
        x_vec = x_vec + gamma0 * eta * w_vec(:, prev)

        res  = abs(sigma0) * res
        eta  = -sigma0 * eta
        beta0 = beta1
        tmp_indx = curr
        curr = prev
        prev = tmp_indx
    enddo

    if (myid == master) then
        write(mystd, "(a25, i5, E10.2, a5, E10.2)") "PMINRES iteration:  ", iter, abs(res), "-->", linsys_tol
    endif

    ! Verify final residual |H x - b|
    tmp_vec = czero
    call pspmv_csr(new_comm, nblock, end_indx, needed, nloc, nloc, ham, x_vec, tmp_vec)
    tmp_vec   = tmp_vec - b_vec
    alpha     = dot_product(tmp_vec, tmp_vec)
    call MPI_ALLREDUCE(alpha, alpha_mpi, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, new_comm, ierror)
    alpha = alpha_mpi
    if (myid == master) then
        write(mystd, "(a34, E10.2)") "Final precision of |H*x-b|:  ", sqrt(abs(alpha))
    endif

    return
end subroutine pminres_csr
