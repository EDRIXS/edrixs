"""Consistency checks for spectra with multiple retained initial states.

These checks exercise complete XAS/RIXS workflows in which several low-energy
states contribute with thermal weights.  They compare the SciPy sparse path
with the existing dense Python spectra.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from edrixs.solvers import rixs, rixs_1v1c_py, xas, xas_1v1c_py

from _helpers import exact_1v1c_reference_data

pytestmark = [
    pytest.mark.integration,
    pytest.mark.filterwarnings(
        "ignore:.*is deprecated; use .* instead.:DeprecationWarning"
    ),
]


def test_xas_two_initial_states_and_energy_dependent_broadening_match_dense(
    small_1v1c_problem,
):
    """Compare full multistate XAS spectra across dense and SciPy workflows.

    The check retains two initial states, applies thermal weighting and several
    polarizations, and uses an incident-energy-dependent lifetime width before
    comparing the final XAS arrays.
    """
    _, hmat_n, transitions, eval_i, evec_i, eval_n, transitions_eig = (
        exact_1v1c_reference_data(small_1v1c_problem)
    )
    kept = [0, 1]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.linspace(center - 0.7, center + 0.7, 7)
    gamma_c = np.linspace(0.18, 0.28, len(ominc))
    pol_type = [("left", 0.0), ("right", 0.0), ("isotropic", 0.0)]

    actual = xas(
        eval_i[kept],
        evec_i[:, kept],
        hmat_n,
        transitions,
        ominc,
        gamma_c=gamma_c,
        thin=0.65,
        phi=0.17,
        pol_type=pol_type,
        temperature=3000.0,
        backend="scipy",
        backend_kws={"nkryl": hmat_n.shape[0]},
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


def test_rixs_two_initial_states_match_dense(small_1v1c_problem):
    """Compare the full multistate RIXS spectrum with the dense legacy solver.

    The check combines two retained states, energy-dependent core and final
    broadenings, and two polarization channels, including the elastic
    contribution in both spectra.
    """
    hmat_i, hmat_n, transitions, eval_i, evec_i, eval_n, transitions_eig = (
        exact_1v1c_reference_data(small_1v1c_problem)
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

    actual = rixs(
        eval_i[kept],
        evec_i[:, kept],
        hmat_i,
        hmat_n,
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
        backend="scipy",
        backend_kws={
            "nkryl": hmat_i.shape[0],
            "linsys_tol": 1e-12,
            "linsys_maxiter": 500,
            "linsys_restart": max(5, hmat_n.shape[0]),
        },
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
        skip_gs=False,
    )

    assert_allclose(actual, expected, rtol=4e-7, atol=4e-9)
