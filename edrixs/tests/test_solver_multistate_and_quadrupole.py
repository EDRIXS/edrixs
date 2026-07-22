import numpy as np
import pytest
from numpy.testing import assert_allclose

import edrixs._solvers_helpers as helpers
from edrixs.solvers import (
    ops,
    rixs_1v1c_py,
    rixs_krylov_scipy,
    xas_1v1c_py,
    xas_krylov_scipy,
)


def _exact_problem_data(problem):
    h_i, h_n, transitions = ops(*problem, backend="dense")
    eval_i, evec_i = np.linalg.eigh(h_i)
    eval_n, evec_n = np.linalg.eigh(h_n)
    transitions_eigenbasis = np.stack(
        [evec_n.conj().T @ transition @ evec_i for transition in transitions]
    )
    return h_i, h_n, transitions, eval_i, evec_i, eval_n, transitions_eigenbasis


@pytest.mark.integration
def test_xas_two_initial_states_and_energy_dependent_broadening_match_dense(
    small_1v1c_problem,
):
    _, h_n, transitions, eval_i, evec_i, eval_n, transitions_eig = (
        _exact_problem_data(small_1v1c_problem)
    )
    kept = [0, 1]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.linspace(center - 0.7, center + 0.7, 7)
    gamma_c = np.linspace(0.18, 0.28, len(ominc))
    pol_type = [("left", 0.0), ("right", 0.0), ("isotropic", 0.0)]

    actual = xas_krylov_scipy(
        eval_i[kept],
        evec_i[:, kept],
        h_n,
        transitions,
        ominc,
        gamma_c=gamma_c,
        thin=0.65,
        phi=0.17,
        pol_type=pol_type,
        temperature=3000.0,
        nkryl=h_n.shape[0],
    )
    expected = xas_1v1c_py(
        eval_i,
        eval_n,
        transitions_eig,
        ominc,
        gamma_c=gamma_c,
        thin=0.65,
        phi=0.17,
        pol_type=pol_type,
        gs_list=kept,
        temperature=3000.0,
    )

    assert_allclose(actual, expected, rtol=3e-8, atol=3e-10)


@pytest.mark.integration
@pytest.mark.parametrize("skip_gs", [False, True])
def test_rixs_two_initial_states_and_skip_gs_match_dense(
    small_1v1c_problem,
    skip_gs,
):
    h_i, h_n, transitions, eval_i, evec_i, eval_n, transitions_eig = (
        _exact_problem_data(small_1v1c_problem)
    )
    kept = [0, 1]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.array([center - 0.15, center + 0.2])
    eloss = np.linspace(-0.25, 1.0, 9)
    gamma_c = np.array([0.22, 0.27])
    gamma_f = np.linspace(0.055, 0.085, len(eloss))
    pol_type = [
        ("left", 0.0, "right", 0.0),
        ("isotropic", 0.0, "linear", 0.2),
    ]

    actual = rixs_krylov_scipy(
        eval_i[kept],
        evec_i[:, kept],
        h_i,
        h_n,
        transitions,
        ominc,
        eloss,
        gamma_c=gamma_c,
        gamma_f=gamma_f,
        thin=0.7,
        thout=1.05,
        phi=0.13,
        pol_type=pol_type,
        temperature=3000.0,
        skip_gs=skip_gs,
        nkryl=h_i.shape[0],
        linsys_tol=1e-12,
        linsys_maxiter=500,
        linsys_restart=max(5, h_n.shape[0]),
        workers=1,
    )
    expected = rixs_1v1c_py(
        eval_i,
        eval_n,
        transitions_eig,
        ominc,
        eloss,
        gamma_c=gamma_c,
        gamma_f=gamma_f,
        thin=0.7,
        thout=1.05,
        phi=0.13,
        pol_type=pol_type,
        gs_list=kept,
        temperature=3000.0,
        skip_gs=skip_gs,
    )

    assert_allclose(actual, expected, rtol=4e-7, atol=4e-9)


def test_skip_gs_projects_every_retained_initial_state(monkeypatch):
    monkeypatch.setattr(
        helpers,
        "_gmres_scipy_compat",
        lambda *args, **kwargs: (np.array([2.0, -3.0], dtype=complex), 0),
    )

    common = dict(
        hmat_i=np.eye(2),
        hmat_n=np.eye(2),
        trans_op_H=[np.eye(2), np.zeros((2, 2)), np.zeros((2, 2))],
        polvec_f=np.array([1.0, 0.0, 0.0]),
        eval_i=np.array([0.0, 0.4]),
        evec_i=np.eye(2),
        istate=0,
        omega=1.0,
        gamma_c=0.2,
        rhs=np.ones(2),
        nkryl=2,
        linsys_tol=1e-10,
        linsys_maxiter=10,
        linsys_restart=5,
    )

    unprojected = helpers._rixs_krylov_one_contribution_scipy(
        **common,
        skip_gs=False,
    )
    projected = helpers._rixs_krylov_one_contribution_scipy(
        **common,
        skip_gs=True,
    )

    assert unprojected["norm"] == pytest.approx(13.0)
    assert projected["norm"] == 0.0
    assert projected["npoles"] == 1
    assert_allclose(projected["alpha"], [0.0])
    assert projected["beta"].size == 0


def _random_hermitian(rng, n):
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (raw + raw.conj().T) / 2


@pytest.mark.integration
def test_quadrupole_xas_matches_dense_reference_on_tiny_matrices():
    rng = np.random.default_rng(41)
    h_n = _random_hermitian(rng, 4)
    eval_n, evec_n = np.linalg.eigh(h_n)
    eval_i = np.array([-0.2, 0.35])
    evec_i = np.eye(2, dtype=complex)
    transitions = [
        rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
        for _ in range(5)
    ]
    transitions_eig = np.stack(
        [evec_n.conj().T @ transition @ evec_i for transition in transitions]
    )
    ominc = np.linspace(-1.0, 1.4, 8)
    pol_type = [("linear", 0.3), ("left", 0.0), ("isotropic", 0.0)]

    actual = xas_krylov_scipy(
        eval_i,
        evec_i,
        h_n,
        transitions,
        ominc,
        gamma_c=0.17,
        thin=0.6,
        phi=0.2,
        pol_type=pol_type,
        temperature=5000.0,
        nkryl=4,
    )
    expected = xas_1v1c_py(
        eval_i,
        eval_n,
        transitions_eig,
        ominc,
        gamma_c=0.17,
        thin=0.6,
        phi=0.2,
        pol_type=pol_type,
        gs_list=[0, 1],
        temperature=5000.0,
    )

    assert_allclose(actual, expected, rtol=5e-9, atol=5e-11)


@pytest.mark.integration
def test_quadrupole_rixs_matches_dense_reference_on_tiny_matrices():
    rng = np.random.default_rng(52)
    h_i = _random_hermitian(rng, 3)
    h_n = _random_hermitian(rng, 4)
    eval_i, evec_i = np.linalg.eigh(h_i)
    eval_n, evec_n = np.linalg.eigh(h_n)
    transitions = [
        rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
        for _ in range(5)
    ]
    transitions_eig = np.stack(
        [evec_n.conj().T @ transition @ evec_i for transition in transitions]
    )
    ominc = np.array([-0.25, 0.35])
    eloss = np.linspace(-0.5, 1.2, 8)
    pol_type = [("linear", 0.2, "left", 0.0)]

    actual = rixs_krylov_scipy(
        eval_i,
        evec_i,
        h_i,
        h_n,
        transitions,
        ominc,
        eloss,
        gamma_c=0.21,
        gamma_f=0.09,
        thin=0.55,
        thout=1.0,
        phi=0.16,
        pol_type=pol_type,
        temperature=5000.0,
        nkryl=3,
        linsys_tol=1e-13,
        linsys_maxiter=200,
        linsys_restart=8,
        workers=1,
    )
    expected = rixs_1v1c_py(
        eval_i,
        eval_n,
        transitions_eig,
        ominc,
        eloss,
        gamma_c=0.21,
        gamma_f=0.09,
        thin=0.55,
        thout=1.0,
        phi=0.16,
        pol_type=pol_type,
        gs_list=[0, 1, 2],
        temperature=5000.0,
    )

    assert_allclose(actual, expected, rtol=2e-8, atol=2e-10)
