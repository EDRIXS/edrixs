import numpy as np
import pytest
import edrixs
from edrixs.soc import atom_hsoc


@pytest.mark.parametrize("case,expected_shape", [
    ('p', (6, 6)),
    ('t2g', (6, 6)),
    ('d', (10, 10)),
    ('f', (14, 14)),
])
def test_soc_shape(case, expected_shape):
    """SOC matrix has correct shape for each shell."""
    hsoc = atom_hsoc(case, 1.0)
    assert hsoc.shape == expected_shape


@pytest.mark.parametrize("case", ['p', 't2g', 'd', 'f'])
def test_soc_hermitian(case):
    """SOC matrix is Hermitian for each shell."""
    hsoc = atom_hsoc(case, 1.0)
    assert np.allclose(hsoc, np.conj(hsoc.T))


@pytest.mark.parametrize("case", ['p', 't2g', 'd', 'f'])
def test_soc_zero_strength_gives_zero_matrix(case):
    """SOC matrix is zero when soc=0."""
    hsoc = atom_hsoc(case, 0.0)
    assert np.allclose(hsoc, 0)


@pytest.mark.parametrize("case", ['p', 't2g', 'd', 'f'])
def test_soc_scales_linearly_with_strength(case):
    """SOC matrix scales linearly with soc parameter."""
    soc1, soc2 = 0.5, 2.0
    hsoc1 = atom_hsoc(case, soc1)
    hsoc2 = atom_hsoc(case, soc2)
    assert np.allclose(hsoc2, hsoc1 * (soc2 / soc1))


@pytest.mark.parametrize("case", ['p', 'd', 'f'])
def test_soc_traceless(case):
    """SOC matrix is traceless (eigenvalues sum to zero)."""
    hsoc = atom_hsoc(case, 1.0)
    assert np.isclose(np.trace(hsoc), 0)


def test_soc_d_eigenvalues():
    """d-shell SOC eigenvalues: 6 states at +soc, 4 states at -3*soc/2."""
    soc = 0.1
    hsoc = atom_hsoc('d', soc)
    evals = np.sort(np.linalg.eigvalsh(hsoc))
    # j=5/2 (6 states): eigenvalue = soc * l·s = soc * 1 = soc
    # j=3/2 (4 states): eigenvalue = soc * l·s = soc * (-3/2) = -3*soc/2
    expected = np.sort([-3 * soc / 2] * 4 + [soc] * 6)
    assert np.allclose(evals, expected, atol=1e-10)


def test_soc_p_eigenvalues():
    """p-shell SOC eigenvalues: 4 states at +soc/2, 2 states at -soc."""
    soc = 0.2
    hsoc = atom_hsoc('p', soc)
    evals = np.sort(np.linalg.eigvalsh(hsoc))
    # For p-shell (l=1, s=1/2):
    # j=3/2: l·s = (J^2-L^2-S^2)/2 = (15/4 - 2 - 3/4)/2 = (8/4)/2 = 1/2
    # j=1/2: l·s = (3/4 - 2 - 3/4)/2 = (-8/4)/2 = -1
    # Eigenvalues: 4 states at soc/2, 2 states at -soc
    expected = np.sort([-soc] * 2 + [soc / 2] * 4)
    assert np.allclose(evals, expected, atol=1e-10)


def test_soc_t2g_eigenvalues():
    """t2g SOC (effective l=1 with negative sign) eigenvalues."""
    soc = 0.15
    hsoc = atom_hsoc('t2g', soc)
    evals = np.sort(np.linalg.eigvalsh(hsoc))
    # t2g uses negative SOC (T-P equivalence: leff=1, opposite sign)
    # j=3/2: 4 states at -soc/2; j=1/2: 2 states at +soc
    expected = np.sort([-soc / 2] * 4 + [soc] * 2)
    assert np.allclose(evals, expected, atol=1e-10)


def test_soc_d_eigenvalue_count():
    """d-shell SOC has 10 eigenvalues (5 orbitals × 2 spins)."""
    hsoc = atom_hsoc('d', 0.1)
    evals = np.linalg.eigvalsh(hsoc)
    assert len(evals) == 10


def test_soc_p_eigenvalue_count():
    """p-shell SOC has 6 eigenvalues (3 orbitals × 2 spins)."""
    hsoc = atom_hsoc('p', 0.1)
    evals = np.linalg.eigvalsh(hsoc)
    assert len(evals) == 6


def test_soc_f_eigenvalue_count():
    """f-shell SOC has 14 eigenvalues (7 orbitals × 2 spins)."""
    hsoc = atom_hsoc('f', 0.1)
    evals = np.linalg.eigvalsh(hsoc)
    assert len(evals) == 14


def test_soc_f_eigenvalues():
    """f-shell SOC eigenvalues: 8 states at +soc*3/2, 6 states at -2*soc."""
    soc = 0.05
    hsoc = atom_hsoc('f', soc)
    evals = np.sort(np.linalg.eigvalsh(hsoc))
    # For f-shell (l=3, s=1/2):
    # j=7/2: l·s = (63/4 - 12 - 3/4)/2 = (36/4)/2 = 9/2... let me recalculate
    # j=7/2: (j(j+1) - l(l+1) - s(s+1))/2 = (63/4 - 12 - 3/4)/2 = (48/4)/2 = 3/2
    # j=5/2: (j(j+1) - l(l+1) - s(s+1))/2 = (35/4 - 12 - 3/4)/2 = (8/4 - 48/4)/2
    #       = (8/4 - 48/4)/2... wait: (35/4 - 48/4 - 3/4)/2 = (-16/4)/2 = -2
    expected = np.sort([-2 * soc] * 6 + [3 * soc / 2] * 8)
    assert np.allclose(evals, expected, atol=1e-10)
