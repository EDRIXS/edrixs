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
@pytest.mark.parametrize('ntrans', [3, 5])
def test_dense_xas_matches_legacy_exact_eigenstate_sum(ntrans):
    """Staged dense XAS agrees with the legacy dense implementation."""
    _, hmat_n, transitions, eval_i, evec_i, eval_n, transitions_eigen = (
        exact_problem(ntrans=ntrans)
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
@pytest.mark.parametrize('ntrans', [3, 5])
@pytest.mark.filterwarnings(
    "ignore:.*is deprecated; use .* instead.:DeprecationWarning"
)
def test_dense_rixs_matches_legacy_exact_eigenstate_sum(skip_gs, ntrans):
    """Staged dense RIXS agrees with the legacy dense implementation."""
    problem = exact_problem(ntrans=ntrans)
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
    kept = [1, 0]
    spectrum, poles = rixs(
        eval_i[kept], evec_i[:, kept], hmat_i, hmat_n, transitions,
        ominc=[0.2, 0.6], eloss=energy_loss, gamma_f=gamma_final,
        pol_type=[('linear', 0, 'linear', 0), ('left', 0, 'right', 0)],
        return_poles=True, temperature=3000.0,
        backend='dense',
    )

    assert spectrum.shape == (2, 2, 2)
    assert len(poles) == 2
    for incident_index, channels in enumerate(poles):
        assert len(channels) == 2
        for channel, record in enumerate(channels):
            assert set(record) == {'npoles', 'eigval', 'norm', 'alpha', 'beta'}
            assert_allclose(record['eigval'], eval_i[kept])
            assert_allclose(
                spectrum[incident_index, :, channel],
                get_spectra_from_poles(
                    record, energy_loss, gamma_final, temperature=3000.0
                ),
                rtol=1e-12,
                atol=1e-12,
            )


@pytest.mark.parametrize('array', [list, np.array])
def test_dense_real_inputs_match_single_level_response(array):
    """Array-like real inputs support both spectra without forced complex casts."""
    energies, vectors = array([0]), array([[1]])
    hi, hn = array([[0]]), array([[2]])
    transitions = array([[[1]], [[1]], [[1]]])
    ominc, eloss = array([1., 2., 3.]), array([-0.1, 0., 0.1])
    denominator = (np.array(ominc) - 2)**2 + 0.2**2

    absorption = xas(
        energies, vectors, hn, transitions, ominc, gamma_c=0.2,
        backend='dense',
    )
    scattering = rixs(
        energies, vectors, hi, hn, transitions, ominc, eloss,
        gamma_c=0.2, gamma_f=0.05, thin=0, thout=0,
        scatter_axis=np.eye(3).tolist(), backend='dense',
    )
    assert_allclose(absorption[:, 0], 0.2 / (np.pi * denominator))
    final_response = 0.05 / (np.pi * (np.array(eloss)**2 + 0.05**2))
    assert_allclose(
        scattering[:, :, 0], final_response[None, :] / denominator[:, None]
    )


@pytest.mark.parametrize('skip_gs', [False, True])
def test_dense_zero_final_vector_returns_zero_poles(skip_gs):
    """Zero transitions and excluded final states produce empty responses."""
    transition = np.ones((1, 1)) if skip_gs else np.zeros((1, 1))
    spectrum, poles = rixs(
        [0], [[1]], [[0]], [[2]], [transition] * 3, [2.0], [0.0],
        thin=0, thout=0, skip_gs=skip_gs, return_poles=True, backend='dense',
    )
    assert_allclose(spectrum, 0)
    assert poles[0][0]['norm'] == [0.0]
    assert poles[0][0]['npoles'] == [1]
    assert_allclose(poles[0][0]['alpha'][0], [0.0])
    assert poles[0][0]['beta'][0].size == 0
