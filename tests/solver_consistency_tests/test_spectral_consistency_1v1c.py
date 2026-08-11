"""Consistency checks for complete small-model 1v1c spectral workflows.

These checks pass a model through operator construction and the public SciPy
XAS/RIXS drivers, then compare complete spectra or outputs with the existing
dense Python path.  They are not isolated unit tests.
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


def test_public_xas_matches_existing_dense_solver(small_1v1c_problem):
    """Compare the complete SciPy and dense 1v1c XAS command chains.

    Both paths start from the same setup output and initial states, apply photon
    transitions, evaluate the intermediate response, apply broadening, and
    return the final incident-energy spectrum for linear and isotropic light.
    """
    _, hmat_n, transitions, eval_i, evec_i, eval_n, trans_eig = (
        exact_1v1c_reference_data(small_1v1c_problem)
    )
    kept = [0]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.linspace(center - 0.8, center + 0.8, 9)
    pol_type = [("linear", 0.2), ("isotropic", 0.0)]

    sparse_result = xas(
        eval_i[kept],
        evec_i[:, kept],
        hmat_n,
        transitions,
        ominc,
        gamma_c=0.25,
        thin=0.7,
        phi=0.1,
        pol_type=pol_type,
        temperature=20.0,
        backend="scipy",
        backend_kws={"nkryl": hmat_n.shape[0]},
    )
    dense_result = xas_1v1c_py(
        eval_i,
        eval_n,
        trans_eig,
        ominc,
        gamma_c=0.25,
        thin=0.7,
        phi=0.1,
        pol_type=pol_type,
        gs_list=kept,
        temperature=20.0,
    )

    assert_allclose(sparse_result, dense_result, rtol=2e-8, atol=2e-10)


def test_public_rixs_matches_existing_dense_solver(small_1v1c_problem):
    """Compare the complete SciPy and dense 1v1c RIXS command chains.

    The check covers incoming transition, intermediate correction-vector solve,
    outgoing transition, final-state response, broadening, and assembly of the
    incident-energy by energy-loss spectrum.
    """
    hmat_i, hmat_n, transitions, eval_i, evec_i, eval_n, trans_eig = (
        exact_1v1c_reference_data(small_1v1c_problem)
    )
    kept = [0]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.array([center - 0.2, center + 0.25])
    eloss = np.linspace(-0.4, 1.2, 11)
    pol_type = [("linear", 0.1, "linear", -0.3)]

    sparse_result = rixs(
        eval_i[kept],
        evec_i[:, kept],
        hmat_i,
        hmat_n,
        transitions,
        ominc,
        eloss,
        gamma_c=0.25,
        gamma_f=0.08,
        thin=0.7,
        thout=1.1,
        phi=0.2,
        pol_type=pol_type,
        temperature=20.0,
        backend="scipy",
        backend_kws={
            "nkryl": hmat_i.shape[0],
            "linsys_tol": 1e-12,
            "linsys_maxiter": 500,
            "linsys_restart": max(5, hmat_n.shape[0]),
        },
    )
    dense_result = rixs_1v1c_py(
        eval_i,
        eval_n,
        trans_eig,
        ominc,
        eloss,
        gamma_c=0.25,
        gamma_f=0.08,
        thin=0.7,
        thout=1.1,
        phi=0.2,
        pol_type=pol_type,
        gs_list=kept,
        temperature=20.0,
    )

    assert_allclose(sparse_result, dense_result, rtol=2e-7, atol=2e-9)


def test_rixs_return_poles_contract(small_1v1c_problem):
    """Check the optional pole output of the complete SciPy RIXS workflow.

    After each correction-vector and final-state response calculation, the
    driver can return the intermediate pole records used to assemble the final
    spectrum.  This check verifies their nesting and fields at the public API.
    """
    hmat_i, hmat_n, transitions, eval_i, evec_i, _, _ = (
        exact_1v1c_reference_data(small_1v1c_problem)
    )
    spectrum, poles = rixs(
        eval_i[:1],
        evec_i[:, :1],
        hmat_i,
        hmat_n,
        transitions,
        ominc=[1.0],
        eloss=[0.0, 0.2],
        gamma_c=0.2,
        gamma_f=0.1,
        return_poles=True,
        backend="scipy",
        backend_kws={
            "nkryl": hmat_i.shape[0],
        },
    )

    assert spectrum.shape == (1, 2, 1)
    assert len(poles) == 1 and len(poles[0]) == 1
    assert set(poles[0][0]) == {"npoles", "eigval", "norm", "alpha", "beta"}
    assert len(poles[0][0]["eigval"]) == 1
