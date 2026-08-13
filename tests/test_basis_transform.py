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


# --- Independent regression oracles ---

def test_tmat_c2r_p_known_coefficients():
    """Pin the p-shell spherical-to-real convention, including phases."""
    s2 = np.sqrt(2.0)
    expected = np.array([
        [1 / s2, 1j / s2, 0],
        [0, 0, 1],
        [-1 / s2, 1j / s2, 0],
    ], dtype=complex)
    assert np.allclose(tmat_c2r('p'), expected)


def test_tmat_c2r_d_known_coefficients():
    """Pin d-orbital ordering and phases, which unitarity alone cannot detect."""
    s2 = np.sqrt(2.0)
    expected = np.zeros((5, 5), dtype=complex)
    expected[2, 0] = 1
    expected[1, 1], expected[3, 1] = 1 / s2, -1 / s2
    expected[1, 2], expected[3, 2] = 1j / s2, 1j / s2
    expected[0, 3], expected[4, 3] = 1 / s2, 1 / s2
    expected[0, 4], expected[4, 4] = 1j / s2, -1j / s2
    assert np.allclose(tmat_c2r('d'), expected)


def test_tmat_r2cub_f_known_coefficients():
    """Pin representative f cubic-harmonic coefficients."""
    T = tmat_r2cub_f()
    assert np.isclose(T[1, 0], -np.sqrt(6.0) / 4)
    assert np.isclose(T[5, 0], np.sqrt(10.0) / 4)
    assert np.isclose(T[0, 2], 1.0)
    assert np.isclose(T[4, 6], 1.0)


def test_tmat_c2j_l1_known_clebsch_gordan_coefficients():
    """Pin representative Clebsch-Gordan coefficients and basis ordering."""
    T = tmat_c2j(1)
    assert np.isclose(T[0, 0], -np.sqrt(2.0 / 3.0))
    assert np.isclose(T[3, 0], np.sqrt(1.0 / 3.0))
    assert np.isclose(T[1, 2], 1.0)
    assert np.isclose(T[4, 5], 1.0)


def test_cb_op_matches_independent_matrix_formula_with_distinct_transforms():
    """Use an explicit oracle rather than another call to cb_op."""
    op = np.array([[1 + 2j, 2 - 1j], [3 + 0.5j, -4j]])
    TL = np.array([[0, 1], [1, 0]], dtype=complex)
    TR = np.array([[1, 0], [0, 1j]], dtype=complex)
    expected = TL.conj().T @ op @ TR
    assert np.allclose(cb_op(op, TL, TR), expected)
    assert np.allclose(cb_op2(op, TL, TR), expected)


def test_cb_op_batch_matches_independent_matrix_formula():
    """Batch behavior is checked against NumPy, not scalar cb_op calls."""
    ops = np.array([
        [[1, 2j], [3, 4]],
        [[-1j, 2], [0.5, 3j]],
    ], dtype=complex)
    T = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    expected = np.stack([T.conj().T @ op @ T for op in ops])
    assert np.allclose(cb_op(ops, T), expected)


def test_transform_utensor_nonidentity_matches_einsum_oracle():
    """A nonidentity rank-4 transform catches a no-op implementation."""
    umat = np.zeros((2, 2, 2, 2), dtype=complex)
    umat[0, 1, 1, 0] = 2.0 + 0.5j
    umat[1, 0, 0, 1] = -1.0j
    T = np.array([[1, 1], [1j, -1j]], dtype=complex) / np.sqrt(2)
    expected = np.einsum('ai,bj,abcd,ck,dl->ijkl',
                         T.conj(), T.conj(), umat, T, T)
    assert np.allclose(transform_utensor(umat, T), expected)
