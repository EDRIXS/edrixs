import numpy as np
import pytest
import edrixs
from edrixs.coulomb_utensor import (
    get_gaunt, umat_slater, get_umat_slater, get_umat_kanamori,
    get_umat_kanamori_ge, get_F0
)


# --- Shape tests ---

@pytest.mark.parametrize("case,expected_n", [
    ('s', 2), ('p', 6), ('t2g', 6), ('d', 10), ('f', 14),
    ('p12', 2), ('p32', 4), ('d32', 4), ('d52', 6),
])
def test_umat_slater_single_shell_shape(case, expected_n):
    """Single-shell Coulomb tensor has shape (n, n, n, n)."""
    if case == 's':
        umat = get_umat_slater(case, 1.0)
    elif case in ('p', 'p12', 'p32'):
        umat = get_umat_slater(case, 1.0, 0.5)
    elif case in ('t2g', 'd', 'd32', 'd52'):
        umat = get_umat_slater(case, 1.0, 0.5, 0.3)
    else:
        umat = get_umat_slater(case, 1.0, 0.5, 0.3, 0.1)
    assert umat.shape == (expected_n, expected_n, expected_n, expected_n)


def test_umat_slater_dp_shape():
    """d+p two-shell Coulomb tensor has shape (16, 16, 16, 16)."""
    umat = get_umat_slater('dp', 4.0, 2.0, 0.5, 0.0, 1.0, 0.5, 0.3, 0.0, 0.0)
    assert umat.shape == (16, 16, 16, 16)


# --- Zero Slater integrals give zero tensor ---

@pytest.mark.parametrize("case,args", [
    ('s', (0.0,)),
    ('p', (0.0, 0.0)),
    ('d', (0.0, 0.0, 0.0)),
    ('t2g', (0.0, 0.0, 0.0)),
])
def test_umat_slater_zero_integrals(case, args):
    """Coulomb tensor is zero when all Slater integrals are zero."""
    umat = get_umat_slater(case, *args)
    assert np.allclose(umat, 0)


# --- Hermitian symmetry: U[i,j,k,l] = conj(U[j,i,l,k]) ---

def test_umat_slater_d_hermitian_symmetry():
    """d-shell Coulomb tensor satisfies U[i,j,k,l] = conj(U[j,i,l,k])."""
    umat = get_umat_slater('d', 4.0, 2.0, 0.5)
    assert np.allclose(umat, np.conj(np.transpose(umat, (1, 0, 3, 2))))


def test_umat_slater_p_hermitian_symmetry():
    """p-shell Coulomb tensor satisfies Hermitian symmetry."""
    umat = get_umat_slater('p', 2.0, 1.0)
    assert np.allclose(umat, np.conj(np.transpose(umat, (1, 0, 3, 2))))


def test_umat_slater_t2g_hermitian_symmetry():
    """t2g Coulomb tensor satisfies Hermitian symmetry."""
    umat = get_umat_slater('t2g', 4.0, 2.0, 0.5)
    assert np.allclose(umat, np.conj(np.transpose(umat, (1, 0, 3, 2))))


# --- get_F0 formulas ---

def test_get_F0_s_is_zero():
    """F0 for s-shell is zero."""
    assert get_F0('s') == 0.0


def test_get_F0_p_formula():
    """F0 for p-shell: F0 = 2/25 * F2."""
    F2 = 3.0
    assert np.isclose(get_F0('p', F2), 2.0 / 25.0 * F2)


def test_get_F0_d_formula():
    """F0 for d-shell: F0 = 2/63 * (F2 + F4)."""
    F2, F4 = 3.0, 1.5
    assert np.isclose(get_F0('d', F2, F4), 2.0 / 63.0 * F2 + 2.0 / 63.0 * F4)


def test_get_F0_f_formula():
    """F0 for f-shell: F0 = 4/195*F2 + 2/143*F4 + 100/5577*F6."""
    F2, F4, F6 = 5.0, 3.0, 2.0
    expected = 4.0 / 195.0 * F2 + 2.0 / 143.0 * F4 + 100.0 / 5577.0 * F6
    assert np.isclose(get_F0('f', F2, F4, F6), expected)


def test_get_F0_dp():
    """F0 for dp cross-shell: F0 = 1/15*G1 + 3/70*G3."""
    G1, G3 = 2.0, 1.0
    assert np.isclose(get_F0('dp', G1, G3), 1.0 / 15.0 * G1 + 3.0 / 70.0 * G3)


def test_get_F0_sp():
    """F0 for sp: F0 = 1/6 * G1."""
    G1 = 1.0
    assert np.isclose(get_F0('sp', G1), 1.0 / 6.0 * G1)


def test_get_F0_sd():
    """F0 for sd: F0 = 1/10 * G2."""
    G2 = 2.0
    assert np.isclose(get_F0('sd', G2), 1.0 / 10.0 * G2)


# --- Kanamori tensor ---

def test_umat_kanamori_shape():
    """Kanamori tensor has correct shape (norbs, norbs, norbs, norbs)."""
    norbs = 6
    umat = get_umat_kanamori(norbs, 4.0, 0.5)
    assert umat.shape == (norbs, norbs, norbs, norbs)


def test_umat_kanamori_on_site_same_orbital_opposite_spin():
    """Kanamori U[0,1,1,0] = U for same orbital, opposite spin."""
    U, J = 4.0, 0.5
    umat = get_umat_kanamori(6, U, J)
    # alpha=0 (band=0,spin=up), beta=1 (band=0,spin=dn)
    assert np.isclose(umat[0, 1, 1, 0], U)


def test_umat_kanamori_inter_orbital_same_spin():
    """Kanamori U[0,2,2,0] = U-3J for different orbital, same spin (density-density)."""
    U, J = 4.0, 0.5
    umat = get_umat_kanamori(6, U, J)
    # alpha=0 (band=0,spin=up), beta=2 (band=1,spin=up): U2 - J = (U-2J) - J = U-3J
    assert np.isclose(umat[0, 2, 2, 0], U - 3 * J)


def test_umat_kanamori_inter_orbital_opposite_spin():
    """Kanamori U[0,3,3,0] = U-2J for different orbital, opposite spin (density-density)."""
    U, J = 4.0, 0.5
    umat = get_umat_kanamori(6, U, J)
    # alpha=0 (band=0,spin=up), beta=3 (band=1,spin=dn): only U2 term contributes
    assert np.isclose(umat[0, 3, 3, 0], U - 2 * J)


def test_umat_kanamori_pair_hopping():
    """Kanamori pair hopping: U[0,1,3,2] = J."""
    U, J = 4.0, 0.5
    umat = get_umat_kanamori(6, U, J)
    # alpha=0 (band=0,up), beta=1 (band=0,dn) → gamma=2 (band=1,up), delta=3 (band=1,dn)
    assert np.isclose(umat[0, 1, 3, 2], J)


def test_umat_kanamori_ge_shape():
    """get_umat_kanamori_ge has correct shape."""
    norbs = 4
    umat = get_umat_kanamori_ge(norbs, 3.0, 2.5, 0.5, 0.5, 0.5)
    assert umat.shape == (norbs, norbs, norbs, norbs)


# --- Gaunt coefficients ---

def test_gaunt_shape():
    """get_gaunt returns array of correct shape."""
    l1, l2 = 2, 1
    g = get_gaunt(l1, l2)
    assert g.shape == (l1 + l2 + 1, 2 * l1 + 1, 2 * l2 + 1)


def test_gaunt_is_real():
    """Gaunt coefficients are real-valued."""
    g = get_gaunt(2, 2)
    assert np.allclose(g.imag, 0)


def test_gaunt_dd_selection_rule():
    """Gaunt C(d,d) is zero for odd k (parity selection rule)."""
    g = get_gaunt(2, 2)
    for k in [1, 3]:
        assert np.allclose(g[k, :, :], 0)


# --- umat_slater low-level ---

def test_umat_slater_low_level_d_shape():
    """Low-level umat_slater for d-shell returns (10, 10, 10, 10) tensor."""
    fk = {(0, 1, 1, 1, 1): 4.0, (2, 1, 1, 1, 1): 2.0, (4, 1, 1, 1, 1): 0.5}
    umat = umat_slater([2], fk)
    assert umat.shape == (10, 10, 10, 10)


def test_umat_slater_low_level_zero_fk():
    """Low-level umat_slater with all-zero fk returns zero tensor."""
    fk = {(0, 1, 1, 1, 1): 0.0, (2, 1, 1, 1, 1): 0.0, (4, 1, 1, 1, 1): 0.0}
    umat = umat_slater([2], fk)
    assert np.allclose(umat, 0)
