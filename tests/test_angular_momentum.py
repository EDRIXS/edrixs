import numpy as np
import pytest
import edrixs


@pytest.mark.parametrize("ll", [0, 1, 2, 3])
def test_lz_diagonal_eigenvalues(ll):
    """Lz is diagonal with eigenvalues -l, -l+1, ..., +l."""
    lz = edrixs.get_lz(ll)
    expected_diag = np.arange(-ll, ll + 1, dtype=float)
    assert np.allclose(np.diag(lz).real, expected_diag)
    off_diag = lz - np.diag(np.diag(lz))
    assert np.allclose(off_diag, 0)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_l2_equals_ll_plus1_identity(ll):
    """L² = Lx² + Ly² + Lz² = l(l+1)·I."""
    lx = edrixs.get_lx(ll)
    ly = edrixs.get_ly(ll)
    lz = edrixs.get_lz(ll)
    L2 = lx @ lx + ly @ ly + lz @ lz
    norbs = 2 * ll + 1
    expected = ll * (ll + 1) * np.eye(norbs, dtype=complex)
    assert np.allclose(L2, expected)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_commutation_lx_ly(ll):
    """[Lx, Ly] = i·Lz."""
    lx, ly, lz = edrixs.get_lx(ll), edrixs.get_ly(ll), edrixs.get_lz(ll)
    assert np.allclose(lx @ ly - ly @ lx, 1j * lz)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_commutation_ly_lz(ll):
    """[Ly, Lz] = i·Lx."""
    lx, ly, lz = edrixs.get_lx(ll), edrixs.get_ly(ll), edrixs.get_lz(ll)
    assert np.allclose(ly @ lz - lz @ ly, 1j * lx)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_commutation_lz_lx(ll):
    """[Lz, Lx] = i·Ly."""
    lx, ly, lz = edrixs.get_lx(ll), edrixs.get_ly(ll), edrixs.get_lz(ll)
    assert np.allclose(lz @ lx - lx @ lz, 1j * ly)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_ladd_lminus_are_adjoint(ll):
    """L⁺ and L⁻ are adjoint: (L⁺)† = L⁻."""
    ladd = edrixs.get_ladd(ll)
    lminus = edrixs.get_lminus(ll)
    assert np.allclose(ladd, np.conj(lminus.T))


@pytest.mark.parametrize("ll", [0, 1, 2, 3])
def test_lx_ly_lz_hermitian(ll):
    """Lx, Ly, Lz are Hermitian matrices."""
    for op in [edrixs.get_lx(ll), edrixs.get_ly(ll), edrixs.get_lz(ll)]:
        assert np.allclose(op, np.conj(op.T))


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_lx_is_half_ladd_plus_lminus(ll):
    """Lx = (L⁺ + L⁻) / 2."""
    ladd = edrixs.get_ladd(ll)
    lminus = edrixs.get_lminus(ll)
    assert np.allclose(edrixs.get_lx(ll), (ladd + lminus) / 2)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_ly_is_minus_i_half_ladd_minus_lminus(ll):
    """Ly = -i(L⁺ - L⁻) / 2."""
    ladd = edrixs.get_ladd(ll)
    lminus = edrixs.get_lminus(ll)
    assert np.allclose(edrixs.get_ly(ll), -1j * (ladd - lminus) / 2)


@pytest.mark.parametrize("ll", [1, 2, 3])
def test_lz_with_spin_diagonal_structure(ll):
    """Lz with ispin=True has each Lz eigenvalue appearing twice (once per spin)."""
    lz_spin = edrixs.get_lz(ll, ispin=True)
    diag = np.diag(lz_spin).real
    expected = np.repeat(np.arange(-ll, ll + 1, dtype=float), 2)
    # Sort both since ordering interleaves up/down
    assert np.allclose(np.sort(diag), np.sort(expected))


@pytest.mark.parametrize("ll", [1, 2])
def test_lx_with_spin_shape(ll):
    """get_lx(ll, ispin=True) has shape (2*(2ll+1), 2*(2ll+1))."""
    lx_spin = edrixs.get_lx(ll, ispin=True)
    norbs = 2 * (2 * ll + 1)
    assert lx_spin.shape == (norbs, norbs)


def test_pauli_matrices_squared_is_identity():
    """Each Pauli matrix squares to the identity."""
    sigma = edrixs.get_pauli()
    I2 = np.eye(2, dtype=complex)
    for i in range(3):
        assert np.allclose(sigma[i] @ sigma[i], I2)


def test_pauli_matrices_traceless():
    """Pauli matrices are traceless."""
    sigma = edrixs.get_pauli()
    for i in range(3):
        assert np.isclose(np.trace(sigma[i]), 0)


def test_pauli_algebra():
    """σx σy = iσz, σy σz = iσx, σz σx = iσy."""
    sigma = edrixs.get_pauli()
    assert np.allclose(sigma[0] @ sigma[1], 1j * sigma[2])
    assert np.allclose(sigma[1] @ sigma[2], 1j * sigma[0])
    assert np.allclose(sigma[2] @ sigma[0], 1j * sigma[1])


def test_pauli_hermitian():
    """Pauli matrices are Hermitian."""
    sigma = edrixs.get_pauli()
    for i in range(3):
        assert np.allclose(sigma[i], np.conj(sigma[i].T))


@pytest.mark.parametrize("ll", [1, 2])
def test_sz_eigenvalues_are_half_integers(ll):
    """Sz has eigenvalues ±1/2, each appearing (2ll+1) times."""
    sz = edrixs.get_sz(ll)
    evals = np.sort(np.linalg.eigvalsh(sz))
    n_orbs = 2 * ll + 1
    expected = np.sort(np.array([-0.5] * n_orbs + [0.5] * n_orbs))
    assert np.allclose(evals, expected)


@pytest.mark.parametrize("ll", [1, 2])
def test_spin_operators_commutation(ll):
    """[Sx, Sy] = i·Sz."""
    sx, sy, sz = edrixs.get_sx(ll), edrixs.get_sy(ll), edrixs.get_sz(ll)
    assert np.allclose(sx @ sy - sy @ sx, 1j * sz)


@pytest.mark.parametrize("ll", [1, 2])
def test_spin_operators_hermitian(ll):
    """Sx, Sy, Sz are Hermitian."""
    for op in [edrixs.get_sx(ll), edrixs.get_sy(ll), edrixs.get_sz(ll)]:
        assert np.allclose(op, np.conj(op.T))


def test_euler_to_rmat_identity_angles():
    """Euler angles (0, 0, 0) produce the identity matrix."""
    rmat = edrixs.euler_to_rmat(0, 0, 0)
    assert np.allclose(rmat, np.eye(3))


def test_euler_to_rmat_is_proper_rotation():
    """Rotation matrix satisfies R·R^T = I and det(R) = +1."""
    rmat = edrixs.euler_to_rmat(0.3, 0.7, 1.1)
    assert np.allclose(rmat @ rmat.T, np.eye(3))
    assert np.isclose(np.linalg.det(rmat), 1.0)


def test_euler_rmat_roundtrip():
    """euler_to_rmat then rmat_to_euler recovers an equivalent rotation."""
    alpha, beta, gamma = 0.4, 0.9, 1.5
    rmat = edrixs.euler_to_rmat(alpha, beta, gamma)
    alpha2, beta2, gamma2 = edrixs.rmat_to_euler(rmat)
    rmat2 = edrixs.euler_to_rmat(alpha2, beta2, gamma2)
    assert np.allclose(rmat, rmat2)


def test_dmat_spinor_unitary():
    """Wigner D-matrix for j=1/2 spinor is unitary."""
    dmat = edrixs.dmat_spinor(0.3, 0.5, 0.7)
    assert np.allclose(dmat @ np.conj(dmat.T), np.eye(2))


def test_dmat_spinor_identity_angles():
    """dmat_spinor at (0, 0, 0) gives identity for zero rotation."""
    dmat = edrixs.dmat_spinor(0.0, 0.0, 0.0)
    assert np.allclose(dmat, np.eye(2))


def test_zx_to_rmat_orthogonal():
    """zx_to_rmat produces a proper orthogonal matrix."""
    rmat = edrixs.zx_to_rmat([0, 0, 1], [1, 0, 0])
    assert np.allclose(rmat @ rmat.T, np.eye(3))
    assert np.isclose(np.linalg.det(rmat), 1.0)


def test_zx_to_rmat_standard_axes():
    """Standard z=[0,0,1], x=[1,0,0] produces identity rotation."""
    rmat = edrixs.zx_to_rmat([0, 0, 1], [1, 0, 0])
    assert np.allclose(rmat, np.eye(3))


def test_cf_cubic_d_eigenvalues():
    """Cubic CF has 6 eigenvalues at -0.4·10Dq and 4 at +0.6·10Dq."""
    ten_dq = 2.0
    cf = edrixs.cf_cubic_d(ten_dq)
    evals = np.sort(np.linalg.eigvalsh(cf))
    expected = np.sort([-0.4 * ten_dq] * 6 + [0.6 * ten_dq] * 4)
    assert np.allclose(evals, expected)


def test_cf_tetragonal_d_reduces_to_cubic_when_distortion_zero():
    """cf_tetragonal_d(ten_dq, 0, 0) equals cf_cubic_d(ten_dq)."""
    ten_dq = 2.0
    cf_tet = edrixs.cf_tetragonal_d(ten_dq, 0.0, 0.0)
    cf_cub = edrixs.cf_cubic_d(ten_dq)
    assert np.allclose(cf_tet, cf_cub)
