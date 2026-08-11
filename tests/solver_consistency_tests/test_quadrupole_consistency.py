"""Consistency checks for five-component quadrupole XAS and RIXS workflows."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from edrixs.solvers import rixs, rixs_1v1c_py, xas, xas_1v1c_py

from ._helpers import random_hermitian

pytestmark = [
    pytest.mark.integration,
    pytest.mark.filterwarnings(
        "ignore:.*is deprecated; use .* instead.:DeprecationWarning"
    ),
]


def test_quadrupole_xas_matches_dense_reference_on_tiny_matrices():
    """Compare the complete five-transition-component XAS command chains.

    The check bypasses physical setup with tiny explicit matrices, then tests
    polarization conversion, transition combination, intermediate response,
    thermal weighting, and final broadening against the dense XAS path.
    """
    rng = np.random.default_rng(41)
    hmat_n = random_hermitian(rng, 4)
    eval_n, evec_n = np.linalg.eigh(hmat_n)
    eval_i = np.array([-0.2, 0.35])
    evec_i = np.eye(2, dtype=complex)
    transitions = [
        rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
        for _ in range(5)
    ]
    transitions_eig = np.stack(
        [
            evec_n.conj().T @ transition @ evec_i
            for transition in transitions
        ]
    )
    ominc = np.linspace(-1.0, 1.4, 8)
    pol_type = [("linear", 0.3), ("left", 0.0), ("isotropic", 0.0)]

    actual = xas(
        eval_i,
        evec_i,
        hmat_n,
        transitions,
        ominc,
        gamma_c=0.17,
        thin=0.6,
        phi=0.2,
        pol_type=pol_type,
        temperature=5000.0,
        backend="scipy",
        backend_kws={"nkryl": 4},
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


def test_quadrupole_rixs_matches_dense_reference_on_tiny_matrices():
    """Compare the complete five-transition-component RIXS command chains.

    The check covers incoming quadrupole polarization, the intermediate solve,
    outgoing quadrupole transition, final-state response, and spectrum assembly
    against the dense Python implementation on exactly diagonalizable matrices.
    """
    rng = np.random.default_rng(52)
    hmat_i = random_hermitian(rng, 3)
    hmat_n = random_hermitian(rng, 4)
    eval_i, evec_i = np.linalg.eigh(hmat_i)
    eval_n, evec_n = np.linalg.eigh(hmat_n)
    transitions = [
        rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
        for _ in range(5)
    ]
    transitions_eig = np.stack(
        [
            evec_n.conj().T @ transition @ evec_i
            for transition in transitions
        ]
    )
    ominc = np.array([-0.25, 0.35])
    eloss = np.linspace(-0.5, 1.2, 8)
    pol_type = [("linear", 0.2, "left", 0.0)]

    actual = rixs(
        eval_i,
        evec_i,
        hmat_i,
        hmat_n,
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
        backend="scipy",
        backend_kws={
            "nkryl": 3,
            "linsys_tol": 1e-13,
            "linsys_maxiter": 200,
            "linsys_restart": 8,
        },
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
