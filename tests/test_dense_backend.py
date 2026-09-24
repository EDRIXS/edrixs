"""Tests for the exact dense staged backend."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from edrixs.poles import get_spectra_from_poles
from edrixs.solvers import (
    ed,
    rixs,
    rixs_1v1c_py,
    xas,
    xas_1v1c_py,
)


def random_hermitian(rng, size):
    """Return a small random Hermitian matrix."""
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return (raw + raw.conj().T) / 2


def exact_problem(ntrans=3):
    """Return Hamiltonians, eigenpairs, and transition operators."""
    rng = np.random.default_rng(731 + ntrans)
    hmat_i = random_hermitian(rng, 3)
    hmat_n = random_hermitian(rng, 4)
    eval_i, evec_i = np.linalg.eigh(hmat_i)
    eval_n, evec_n = np.linalg.eigh(hmat_n)
    transitions = np.asarray([
        rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
        for _ in range(ntrans)
    ])
    transitions_eigen = np.einsum(
        'an,kab,bi->kni',
        evec_n.conj(), transitions, evec_i, optimize=True,
    )
    return (
        hmat_i, hmat_n, transitions,
        eval_i, evec_i, eval_n, transitions_eigen,
    )


def test_dense_ed_uses_exact_hermitian_diagonalization():
    """Dense ED returns the requested lowest exact eigenpairs."""
    rng = np.random.default_rng(91)
    hamiltonian = random_hermitian(rng, 6)

    eigenvalues, eigenvectors = ed(
        hamiltonian, num_evals=3, backend='dense'
    )

    assert_allclose(eigenvalues, np.linalg.eigvalsh(hamiltonian)[:3])
    assert_allclose(
        hamiltonian @ eigenvectors,
        eigenvectors * eigenvalues[None, :],
        atol=1e-12,
    )


@pytest.mark.filterwarnings(
    "ignore:.*is deprecated; use .* instead.:DeprecationWarning"
)
def test_dense_xas_matches_legacy_exact_eigenstate_sum():
    """Staged dense XAS agrees with the legacy dense implementation."""
    _, hmat_n, transitions, eval_i, evec_i, eval_n, transitions_eigen = (
        exact_problem(ntrans=5)
    )
    kept = [0, 1]
    ominc = np.linspace(-0.8, 1.3, 9)
    gamma_c = np.linspace(0.14, 0.23, len(ominc))
    pol_type = [('linear', 0.24), ('left', 0.0), ('isotropic', 0.0)]

    actual = xas(
        eval_i[kept], evec_i[:, kept], hmat_n, transitions, ominc,
        gamma_c=gamma_c, thin=0.57, phi=0.19,
        pol_type=pol_type, temperature=4200.0, backend='dense',
    )
    expected = xas_1v1c_py(
        eval_i, eval_n, transitions_eigen, ominc,
        gamma_c=gamma_c, thin=0.57, phi=0.19,
        pol_type=pol_type, gs_list=kept, temperature=4200.0,
    )

    assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)


@pytest.mark.parametrize('skip_gs', [False, True])
@pytest.mark.filterwarnings(
    "ignore:.*is deprecated; use .* instead.:DeprecationWarning"
)
def test_dense_rixs_matches_legacy_exact_eigenstate_sum(skip_gs):
    """Staged dense RIXS agrees with the legacy dense implementation."""
    problem = exact_problem(ntrans=3)
    hmat_i, hmat_n, transitions = problem[:3]
    eval_i, evec_i, eval_n, transitions_eigen = problem[3:]
    kept = [0, 1]
    ominc = np.array([-0.31, 0.27])
    eloss = np.linspace(-0.4, 1.2, 10)
    gamma_c = np.array([0.18, 0.25])
    gamma_f = np.linspace(0.05, 0.09, len(eloss))
    pol_type = [
        ('left', 0.0, 'right', 0.0),
        ('linear', -0.21, 'linear', 0.17),
    ]

    actual = rixs(
        eval_i[kept], evec_i[:, kept], hmat_i, hmat_n, transitions,
        ominc, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        thin=0.63, thout=1.07, phi=0.11, pol_type=pol_type,
        temperature=3500.0, skip_gs=skip_gs, backend='dense',
    )
    expected = rixs_1v1c_py(
        eval_i, eval_n, transitions_eigen, ominc, eloss,
        gamma_c=gamma_c, gamma_f=gamma_f,
        thin=0.63, thout=1.07, phi=0.11, pol_type=pol_type,
        gs_list=kept, temperature=3500.0, skip_gs=skip_gs,
    )

    assert_allclose(actual, expected, rtol=3e-13, atol=3e-13)


def test_dense_rixs_preserves_optional_pole_contract():
    """Exact dense RIXS can return full-dimension pole records."""
    hmat_i, hmat_n, transitions, eval_i, evec_i, _, _ = exact_problem()

    energy_loss = np.array([0.0, 0.3])
    gamma_final = np.array([0.08, 0.11])
    spectrum, poles = rixs(
        eval_i[:1], evec_i[:, :1], hmat_i, hmat_n, transitions,
        ominc=[0.2], eloss=energy_loss, gamma_f=gamma_final,
        return_poles=True,
        backend='dense',
    )

    assert spectrum.shape == (1, 2, 1)
    assert len(poles) == 1 and len(poles[0]) == 1
    assert set(poles[0][0]) == {'npoles', 'eigval', 'norm', 'alpha', 'beta'}
    assert len(poles[0][0]['eigval']) == 1
    assert_allclose(
        spectrum[0, :, 0],
        get_spectra_from_poles(
            poles[0][0], energy_loss, gamma_final, temperature=1.0
        ),
        rtol=1e-12,
        atol=1e-12,
    )
