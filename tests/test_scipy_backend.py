"""Unit tests for private numerical handoffs inside the SciPy backend."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse.linalg import aslinearoperator

import edrixs.scipy_backend as backend

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

def test_zero_rhs_short_circuits_rixs_contribution(monkeypatch):
    """A zero incoming-transition vector must bypass the GMRES solve."""

    def forbidden_gmres(*args, **kwargs):
        raise AssertionError("GMRES must not run for a zero RHS")

    monkeypatch.setattr(backend, "_gmres_scipy_compat", forbidden_gmres)
    record = backend._rixs_krylov_one_contribution_scipy(
        hmat_i=aslinearoperator(np.eye(2)),
        hmat_n=aslinearoperator(np.eye(3)),
        trans_op_H=[],
        polvec_f=np.zeros(0),
        eval_i=np.array([0.0]),
        istate=0,
        omega=1.0,
        gamma_c=0.2,
        rhs=np.zeros(3),
        nkryl=2,
        linsys_tol=1e-10,
        linsys_maxiter=10,
        linsys_restart=5,
    )

    assert record["norm"] == 0
    assert record["npoles"] == 1
    assert_allclose(record["alpha"], [0.0])
    assert record["beta"].size == 0

def test_gmres_failure_is_reported_with_context(monkeypatch):
    """A failed correction-vector solve identifies its RIXS contribution."""
    monkeypatch.setattr(
        backend,
        "_gmres_scipy_compat",
        lambda *args, **kwargs: (np.zeros(2, dtype=complex), 7),
    )

    with pytest.raises(RuntimeError, match=r"istate=0.*omega=1.25.*info=7"):
        backend._rixs_krylov_one_contribution_scipy(
            hmat_i=aslinearoperator(np.eye(2)),
            hmat_n=aslinearoperator(np.eye(2)),
            trans_op_H=[aslinearoperator(np.eye(2))] * 3,
            polvec_f=np.array([1.0, 0.0, 0.0]),
            eval_i=np.array([0.0]),
            istate=0,
            omega=1.25,
            gamma_c=0.2,
            rhs=np.ones(2),
            nkryl=2,
            linsys_tol=1e-10,
            linsys_maxiter=10,
            linsys_restart=5,
        )

def test_rixs_backend_forwards_to_krylov_implementation(monkeypatch):
    """The SciPy RIXS backend forwards numerical options to the Krylov solver."""
    sentinel = object()
    received = {}

    def fake_krylov(*args, **kwargs):
        received.update(kwargs)
        return sentinel

    monkeypatch.setattr(backend, "rixs_krylov_scipy", fake_krylov)

    result = backend.rixs_scipy(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        np.eye(2),
        np.eye(2),
        [np.eye(2)] * 3,
        np.array([1.0]),
        np.array([0.0]),
    )

    assert result is sentinel
    assert received["return_poles"] is False
