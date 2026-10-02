"""Unit tests for SciPy spectra and numerical solver handoffs."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose
from scipy.sparse.linalg import aslinearoperator

import edrixs.scipy_backend.scipy_backend as backend


@pytest.mark.parametrize('array', [list, np.array])
@pytest.mark.parametrize('operator', [np.array, sp.csr_matrix, aslinearoperator])
def test_real_spectroscopy_inputs_match_single_level_response(array, operator):
    """Lists and real arrays retain the analytic XAS and RIXS response."""
    eval_i, evec_i = array([0]), array([[1]])
    hmat_i = operator(np.array([[0]]))
    hmat_n = operator(np.array([[2]]))
    trans_op = [operator(np.array([[1]])) for _ in range(3)]
    ominc = array([1.0, 2.0, 3.0])
    eloss = array([-0.1, 0.0, 0.1])
    gamma_c, gamma_f = 0.2, 0.05
    denominator = (np.array(ominc) - 2)**2 + gamma_c**2

    xas = backend.xas_scipy(
        eval_i, evec_i, hmat_n, trans_op, ominc, gamma_c=gamma_c,
    )
    rixs = backend.rixs_scipy(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
        gamma_c=gamma_c, gamma_f=gamma_f, thin=0, thout=0,
        scatter_axis=np.eye(3).tolist(),
    )

    assert_allclose(xas[:, 0], gamma_c / (np.pi * denominator))
    final_response = gamma_f / (np.pi * (np.array(eloss)**2 + gamma_f**2))
    assert_allclose(rixs[:, :, 0], final_response[None, :] / denominator[:, None])


@pytest.mark.parametrize('operator', [np.array, sp.csr_matrix, sp.csc_array])
@pytest.mark.parametrize('offset_i,offset_n', [
    (2.**40, 2.**40), (-2.**40, -2.**40),
    (-6 * 2.**30, -5 * 2.**30), (0., 2.**40), (2.**40, 0.),
])
def test_spectra_are_invariant_to_large_energy_offsets(operator, offset_i, offset_n):
    """Resonant centering preserves small features despite large core energies."""
    hi = np.diag([0., .25, .75, 1.5])
    hn = np.array([
        [2., .125, .0625, 0], [.125, 2.5, .25, .125],
        [.0625, .25, 3., .125], [0, .125, .125, 3.5],
    ])
    rng = np.random.default_rng(41)
    transitions = [
        (rng.integers(-8, 9, (4, 4)) + 1j * rng.integers(-8, 9, (4, 4))) / 8
        for _ in range(3)
    ]
    energies, vectors = np.array([.25, 0.]), np.eye(4)[:, [1, 0]]
    # Dyadic entries stay exact when offsets are added to the input data.
    # Non-dyadic photons in the common-offset case also test subtraction order.
    ominc = np.array([2.2, 2.7] if offset_i == offset_n else [2.25, 2.75])
    eloss = np.linspace(-.1, 1.7, 31)
    gamma_c = [.125, .25]
    gamma_f = np.linspace(.04, .06, len(eloss))
    hi_shifted = operator(hi + offset_i * np.eye(4))
    hn_shifted = operator(hn + offset_n * np.eye(4))
    ominc_shifted = ominc + (offset_n - offset_i)

    xas_expected = backend.xas_scipy(
        energies, vectors, operator(hn), transitions, ominc,
        gamma_c=gamma_c, temperature=3000., backend_kws={'nkryl': 4},
    )
    xas_actual = backend.xas_scipy(
        energies + offset_i, vectors, hn_shifted, transitions, ominc_shifted,
        gamma_c=gamma_c, temperature=3000., backend_kws={'nkryl': 4},
    )
    assert_allclose(xas_actual, xas_expected, rtol=2e-11, atol=2e-12)

    options = dict(
        gamma_c=gamma_c, gamma_f=gamma_f, temperature=3000.,
        backend_kws={'nkryl': 4, 'linsys_tol': 1e-11,
                     'linsys_maxiter': 20, 'linsys_restart': 4},
    )
    rixs_expected = backend.rixs_scipy(
        energies, vectors, operator(hi), operator(hn), transitions,
        ominc, eloss, **options,
    )
    rixs_actual, poles = backend.rixs_scipy(
        energies + offset_i, vectors, hi_shifted, hn_shifted, transitions,
        ominc_shifted, eloss, return_poles=True, **options,
    )
    assert_allclose(rixs_actual, rixs_expected, rtol=2e-11, atol=2e-12)
    for channels in poles:
        assert_allclose(channels[0]['eigval'], energies + offset_i, rtol=0, atol=0)
        for alpha in channels[0]['alpha']:
            assert np.all(np.abs(alpha - offset_i) < 2)
    assert_allclose(hi_shifted @ np.eye(4), hi + offset_i * np.eye(4), rtol=0, atol=0)
    assert_allclose(hn_shifted @ np.eye(4), hn + offset_n * np.eye(4), rtol=0, atol=0)


def test_apply_linear_combination_handles_complex_and_zero_coefficients():
    """Polarization-weighted operator application supports complex weights."""
    operators = [
        aslinearoperator(np.eye(2)),
        aslinearoperator(np.diag([2.0, 3.0])),
    ]
    vector = np.array([1.0, -2.0])

    actual = backend._apply_linear_combination(
        operators, [1j, -0.5], vector
    )
    expected = 1j * vector - 0.5 * (np.diag([2.0, 3.0]) @ vector)

    assert_allclose(actual, expected)
    assert_allclose(
        backend._apply_linear_combination(operators, [0, 0], vector),
        0,
    )


def test_rixs_poles_preserve_energy_polarization_and_state_order():
    """Returned poles belong to their spectrum slice and input states."""
    energies = np.array([0.4, 0.0])
    ominc, eloss = [1.0, 2.5], [-0.1, 0.0, 0.2]
    gamma_c, gamma_f = [0.2, 0.3], np.array([0.05, 0.06, 0.07])
    temperature = 3000.0
    spectrum, poles = backend.rixs_scipy(
        energies, np.eye(2)[:, ::-1], np.diag([0.0, 0.4]), np.diag([2.0, 3.0]),
        [np.diag([1.0, 2.0]), np.diag([3.0, 4.0]), np.diag([5.0, 6.0])],
        ominc, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        thin=0, thout=0, temperature=temperature, return_poles=True,
        pol_type=[('linear', 0, 'linear', 0),
                  ('linear', np.pi / 2, 'linear', np.pi / 2)],
    )

    assert spectrum.shape == (2, 3, 2)
    assert len(poles) == 2
    for incident_index, omega in enumerate(ominc):
        assert len(poles[incident_index]) == 2
        for polarization_index, amplitudes in enumerate(([6., 5.], [4., 3.])):
            record = poles[incident_index][polarization_index]
            denominator = (
                (omega + energies - [3.0, 2.0])**2 + gamma_c[incident_index]**2
            )
            assert_allclose(record['eigval'], energies)
            assert_allclose(
                record['norm'], np.array(amplitudes)**4 / denominator
            )
            assert_allclose([alpha[0] for alpha in record['alpha']], energies)
            for beta in record['beta']:
                assert_allclose(beta, 0, atol=1e-14)
            assert_allclose(
                spectrum[incident_index, :, polarization_index],
                backend.get_spectra_from_poles(
                    record, eloss, gamma_f, temperature
                ),
            )


def test_zero_rhs_short_circuits_rixs_contribution(monkeypatch):
    """A zero incoming-transition vector must bypass the GMRES solve."""

    def forbidden_gmres(*args, **kwargs):
        raise AssertionError("GMRES must not run for a zero RHS")

    monkeypatch.setattr(backend, "gmres", forbidden_gmres)
    spectrum, poles = backend.rixs_scipy(
        [0.0], [[1.0], [0.0]], np.eye(2), np.eye(3),
        [np.zeros((3, 2))] * 3, [1.0], [0.0, 0.2], return_poles=True,
    )

    assert_allclose(spectrum, 0)
    assert poles[0][0]["norm"] == [0.0]
    assert poles[0][0]["npoles"] == [1]
    assert_allclose(poles[0][0]["alpha"][0], [0.0])
    assert poles[0][0]["beta"][0].size == 0


def test_gmres_failure_is_reported_with_context(monkeypatch):
    """A failed correction-vector solve identifies its RIXS contribution."""
    monkeypatch.setattr(
        backend,
        "gmres",
        lambda *args, **kwargs: (np.zeros(2, dtype=complex), 7),
    )

    with pytest.raises(RuntimeError, match=r"istate=0.*omega=1.25.*info=7"):
        backend.rixs_scipy(
            [0.0], [[1.0], [0.0]], np.eye(2), np.eye(2),
            [np.eye(2)] * 3, [1.25], [0.0],
        )


@pytest.mark.parametrize('tolerance_name', ['tol', 'rtol'])
def test_rixs_backend_uses_gmres_options(monkeypatch, tolerance_name):
    """Both SciPy tolerance conventions receive the requested options."""
    received = []

    def old_gmres(operator, rhs, *, tol, restart, maxiter):
        received.append((tol, restart, maxiter))
        return np.linalg.solve(operator @ np.eye(len(rhs)), rhs), 0

    def new_gmres(operator, rhs, *, rtol, atol, restart, maxiter):
        assert atol == 0.0
        received.append((rtol, restart, maxiter))
        return np.linalg.solve(operator @ np.eye(len(rhs)), rhs), 0

    monkeypatch.setattr(
        backend, "gmres", old_gmres if tolerance_name == 'tol' else new_gmres
    )

    result = backend.rixs_scipy(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        np.eye(2),
        np.eye(2),
        [np.eye(2)] * 3,
        np.array([1.0]),
        np.array([0.0]),
        backend_kws={'linsys_tol': 1e-8, 'linsys_maxiter': 20,
                     'linsys_restart': 2, 'nkryl': 2},
    )

    assert isinstance(result, np.ndarray)
    assert np.all(np.isfinite(result))
    assert received
    assert all(options == (1e-8, 2, 20) for options in received)
