"""Tests for the backend-neutral public solver interface."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import ed, get_H, rixs, xas


def test_get_h_rejects_unknown_backend():
    """The public operator constructor validates backend names."""
    with pytest.raises(ValueError, match="Unknown backend"):
        get_H(np.eye(2), np.zeros((2, 2, 2, 2)), object(), backend="bad")


def test_ed_rejects_nonmapping_backend_options():
    """Public wrappers require ``backend_kws`` to be a mapping or ``None``."""
    with pytest.raises(TypeError, match="mapping"):
        ed(np.eye(2), backend="scipy", backend_kws=[("tol", 1e-8)])


def test_public_xas_uses_inferred_scipy_backend():
    """XAS infers SciPy ownership from all supplied operators."""
    hmat_n = sp.diags([1.0, 2.0], format="csr")
    transitions = [
        sp.eye(2, format="csr"),
        sp.csr_matrix((2, 2)),
        sp.csr_matrix((2, 2)),
    ]

    spectrum = xas(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        hmat_n,
        transitions,
        np.array([0.5, 1.0]),
        pol_type=[("isotropic", 0.0)],
        backend_kws={"nkryl": 2},
    )

    assert spectrum.shape == (2, 1)
    assert np.all(np.isfinite(spectrum))


def test_public_rixs_parallel_workers_one_matches_serial():
    """Parallel selection through ``backend_kws`` preserves serial results."""
    hmat_i = sp.diags([0.0, 0.8], format="csr")
    hmat_n = sp.diags([1.5, 2.2], format="csr")
    transitions = [
        sp.eye(2, format="csr"),
        sp.csr_matrix((2, 2)),
        sp.csr_matrix((2, 2)),
    ]
    common = {
        "gamma_c": 0.2,
        "gamma_f": 0.1,
        "pol_type": [("linear", 0.0, "linear", 0.0)],
        "backend": "scipy",
    }
    solver_options = {
        "nkryl": 2,
        "linsys_tol": 1e-10,
        "linsys_maxiter": 20,
        "linsys_restart": 2,
    }

    serial = rixs(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        hmat_i,
        hmat_n,
        transitions,
        np.array([1.0]),
        np.array([0.0, 0.5]),
        backend_kws=solver_options,
        **common,
    )
    parallel = rixs(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        hmat_i,
        hmat_n,
        transitions,
        np.array([1.0]),
        np.array([0.0, 0.5]),
        backend_kws={**solver_options, "parallel": True, "workers": 1},
        **common,
    )

    assert_allclose(parallel, serial)
