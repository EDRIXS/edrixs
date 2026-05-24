import numpy as np
import pytest
from edrixs.basis_transform import (
    cb_op, cb_op2, tmat_c2r, tmat_r2c, tmat_c2j,
    tmat_r2cub_f, tmat_cub2r_f, transform_utensor
)


@pytest.mark.parametrize("case", ['p', 't2g', 'd', 'f'])
def test_tmat_c2r_unitary(case):
    """tmat_c2r produces a unitary matrix: T†T = I."""
    T = tmat_c2r(case)
    n = T.shape[0]
    assert np.allclose(np.conj(T.T) @ T, np.eye(n))


@pytest.mark.parametrize("case", ['p', 't2g', 'd', 'f'])
def test_tmat_c2r_and_r2c_are_inverses(case):
    """tmat_r2c is the inverse of tmat_c2r: Tc2r @ Tr2c = I."""
    Tc2r = tmat_c2r(case)
    Tr2c = tmat_r2c(case)
    n = Tc2r.shape[0]
    assert np.allclose(Tc2r @ Tr2c, np.eye(n))
    assert np.allclose(Tr2c @ Tc2r, np.eye(n))


@pytest.mark.parametrize("case", ['p', 'd', 'f'])
def test_tmat_c2r_with_spin_unitary(case):
    """tmat_c2r with ispin=True is unitary."""
    T = tmat_c2r(case, ispin=True)
    n = T.shape[0]
    assert np.allclose(np.conj(T.T) @ T, np.eye(n))


@pytest.mark.parametrize("case", ['p', 'd', 'f'])
def test_tmat_r2c_with_spin_is_conjugate_transpose(case):
    """tmat_r2c(ispin=True) equals conjugate transpose of tmat_c2r(ispin=True)."""
    Tc2r = tmat_c2r(case, ispin=True)
    Tr2c = tmat_r2c(case, ispin=True)
    assert np.allclose(Tr2c, np.conj(Tc2r.T))


def test_tmat_r2cub_f_unitary():
    """tmat_r2cub_f is a unitary matrix."""
    T = tmat_r2cub_f()
    assert T.shape == (7, 7)
    assert np.allclose(np.conj(T.T) @ T, np.eye(7))


def test_tmat_cub2r_f_is_inverse_of_r2cub():
    """tmat_cub2r_f is the inverse of tmat_r2cub_f."""
    T1 = tmat_r2cub_f()
    T2 = tmat_cub2r_f()
    assert np.allclose(T1 @ T2, np.eye(7))
    assert np.allclose(T2 @ T1, np.eye(7))


@pytest.mark.parametrize("orb_l", [1, 2, 3])
def test_tmat_c2j_unitary(orb_l):
    """tmat_c2j is unitary for l=1, 2, 3."""
    T = tmat_c2j(orb_l)
    n = T.shape[0]
    assert np.allclose(np.conj(T.T) @ T, np.eye(n))


@pytest.mark.parametrize("orb_l,expected_dim", [(1, 6), (2, 10), (3, 14)])
def test_tmat_c2j_shape(orb_l, expected_dim):
    """tmat_c2j has correct shape for each orbital quantum number."""
    T = tmat_c2j(orb_l)
    assert T.shape == (expected_dim, expected_dim)


def test_cb_op_identity_leaves_operator_unchanged():
    """cb_op with identity transformation leaves operator unchanged."""
    n = 5
    np.random.seed(42)
    A = np.random.random((n, n)) + 1j * np.random.random((n, n))
    op = A + np.conj(A.T)
    op_transformed = cb_op(op, np.eye(n, dtype=complex))
    assert np.allclose(op_transformed, op)


def test_cb_op_preserves_eigenvalues():
    """cb_op with a unitary transformation preserves eigenvalues."""
    n = 5
    np.random.seed(42)
    A = np.random.random((n, n)) + 1j * np.random.random((n, n))
    op = A + np.conj(A.T)
    Q, _ = np.linalg.qr(np.random.random((n, n)) + 1j * np.random.random((n, n)))
    op_transformed = cb_op(op, Q)
    evals_orig = np.sort(np.linalg.eigvalsh(op))
    evals_new = np.sort(np.linalg.eigvalsh(op_transformed))
    assert np.allclose(evals_orig, evals_new)


def test_cb_op_batch_applies_per_matrix():
    """cb_op on a batch of matrices applies the transformation to each."""
    n = 4
    np.random.seed(0)
    Q, _ = np.linalg.qr(np.random.random((n, n)) + 1j * np.random.random((n, n)))
    op_batch = np.random.random((3, n, n)) + 1j * np.random.random((3, n, n))
    result = cb_op(op_batch, Q)
    assert result.shape == (3, n, n)
    for i in range(3):
        expected = cb_op(op_batch[i], Q)
        assert np.allclose(result[i], expected)


def test_cb_op2_equivalent_to_cb_op_when_TR_equals_TL():
    """cb_op2 with TR=TL produces same result as cb_op."""
    n = 4
    np.random.seed(1)
    A = np.random.random((n, n)) + 1j * np.random.random((n, n))
    op = A + np.conj(A.T)
    Q, _ = np.linalg.qr(np.random.random((n, n)) + 1j * np.random.random((n, n)))
    assert np.allclose(cb_op(op, Q), cb_op2(op, Q, Q))


def test_transform_utensor_identity_leaves_tensor_unchanged():
    """transform_utensor with identity matrix gives same tensor."""
    n = 4
    umat = np.zeros((n, n, n, n), dtype=complex)
    umat[0, 1, 1, 0] = 1.0
    umat[1, 0, 0, 1] = 1.0
    umat_new = transform_utensor(umat, np.eye(n, dtype=complex))
    assert np.allclose(umat_new, umat)


def test_tmat_c2r_d_has_correct_real_basis_structure():
    """tmat_c2r for d-shell maps complex harmonics to real d-orbitals correctly."""
    T = tmat_c2r('d')
    # T maps from complex to real: each real orbital is a superposition
    # of complex ones. The matrix should be 5×5.
    assert T.shape == (5, 5)
    # The transformation must be unitary
    assert np.allclose(np.conj(T.T) @ T, np.eye(5))


def test_tmat_c2r_round_trip_operator():
    """Applying c2r then r2c basis change is the identity operation."""
    ll = 2
    norbs = 2 * ll + 1
    np.random.seed(5)
    A = np.random.random((norbs, norbs)) + 1j * np.random.random((norbs, norbs))
    op = A + np.conj(A.T)
    Tc2r = tmat_c2r('d')
    Tr2c = tmat_r2c('d')
    op_prime = cb_op(cb_op(op, Tc2r), Tr2c)
    assert np.allclose(op_prime, op)
