"""Tests for the backend-neutral public solver interface."""

import numpy as np
import pytest
import scipy.sparse as sp

from edrixs.solvers import build_op, ed, rixs, xas


def test_build_op_rejects_unknown_backend():
    """The public operator constructor validates backend names."""
    with pytest.raises(ValueError, match="Unknown backend"):
        build_op(np.eye(2), np.zeros((2, 2, 2, 2)), object(), backend="bad")


@pytest.mark.parametrize("backend", [" SciPy ", "SCIPY", "DENSE"])
def test_build_op_requires_exact_backend_name(backend):
    """Backend selectors use the exact documented lowercase names."""
    with pytest.raises(ValueError, match="Unknown backend"):
        build_op(np.eye(2), np.zeros((2, 2, 2, 2)), object(), backend=backend)


def test_ed_rejects_nonmapping_backend_options():
    """Public wrappers require ``backend_kws`` to be a mapping or ``None``."""
    with pytest.raises(TypeError, match="mapping"):
        ed(np.eye(2), backend="scipy", backend_kws=[("tol", 1e-8)])


def test_solver_entry_points_reject_unknown_explicit_backend():
    """Every public solver rejects selectors outside the documented names."""
    with pytest.raises(ValueError, match="Unknown backend"):
        ed(np.eye(2), backend="SCIPY")
    with pytest.raises(ValueError, match="Unknown backend"):
        xas(None, None, None, [], None, backend="SCIPY")
    with pytest.raises(ValueError, match="Unknown backend"):
        rixs(None, None, None, None, [], None, None, backend="SCIPY")


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


def test_public_rixs_uses_inferred_scipy_backend():
    """RIXS infers the SciPy backend from the supplied operators."""
    hmat_i = sp.diags([0.0, 0.8], format="csr")
    hmat_n = sp.diags([1.5, 2.2], format="csr")
    transitions = [
        sp.eye(2, format="csr"),
        sp.csr_matrix((2, 2)),
        sp.csr_matrix((2, 2)),
    ]

    spectrum = rixs(
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        hmat_i,
        hmat_n,
        transitions,
        np.array([1.0]),
        np.array([0.0, 0.5]),
        gamma_c=0.2,
        gamma_f=0.1,
        pol_type=[("linear", 0.0, "linear", 0.0)],
        backend_kws={
            "nkryl": 2,
            "linsys_tol": 1e-10,
            "linsys_maxiter": 20,
            "linsys_restart": 2,
        },
    )

    assert spectrum.shape == (1, 2, 1)
    assert np.all(np.isfinite(spectrum))


def test_public_rixs_skip_gs_removes_retained_final_states():
    """skip_gs removes transitions into the retained initial-state subspace."""
    hmat_i = sp.diags([0.0, 0.8], format="csr")
    hmat_n = sp.diags([1.5, 2.2], format="csr")
    transitions = [sp.eye(2, format="csr")] * 3
    arguments = (
        np.array([0.0]),
        np.array([[1.0], [0.0]]),
        hmat_i,
        hmat_n,
        transitions,
        np.array([1.0]),
        np.array([0.0, 0.5]),
    )
    keywords = {
        "gamma_c": 0.2,
        "gamma_f": 0.1,
        "pol_type": [("linear", np.pi / 2, "linear", np.pi / 2)],
        "backend_kws": {
            "nkryl": 2,
            "linsys_tol": 1e-10,
            "linsys_maxiter": 20,
            "linsys_restart": 2,
        },
    }

    included = rixs(*arguments, **keywords)
    skipped = rixs(*arguments, skip_gs=True, **keywords)

    assert np.any(included > 0)
    np.testing.assert_allclose(skipped, 0.0, atol=1e-14)
