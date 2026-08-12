"""Unit tests for backend-neutral helpers used by :mod:`edrixs.solvers`."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

import edrixs._solvers_helpers as helpers


def test_dense_to_sparse_u_preserves_flattening_convention():
    """Dense and flattened sparse Coulomb data must use the same ordering."""
    rng = np.random.default_rng(3)
    umat = rng.normal(size=(3, 3, 3, 3))
    umat[np.abs(umat) < 0.6] = 0

    actual = helpers._umat_dense_to_sparse(umat).toarray()

    assert_allclose(actual, umat.reshape(9, 9))


def test_sparse_siam_embedding_matches_dense_embedding():
    """Sparse SIAM embedding must match the dense reference placement."""
    rng = np.random.default_rng(8)
    compact = rng.normal(size=(4, 4, 4, 4))
    dense = helpers._embed_impurity_core_umat(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )
    sparse = helpers._embed_impurity_core_umat_sparse(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )

    assert_allclose(sparse.toarray(), dense.reshape(64, 64))


def test_expand_broadening_accepts_scalar_or_exact_length_array():
    """Broadening normalization accepts a scalar or one value per grid point."""
    assert_allclose(helpers._expand_broadening(0.2, 3, "gamma"), [0.2] * 3)
    assert_allclose(
        helpers._expand_broadening([0.1, 0.2], 2, "gamma"),
        [0.1, 0.2],
    )

    with pytest.raises(ValueError, match="shape"):
        helpers._expand_broadening([0.1], 2, "gamma")


def test_infer_backend_recognizes_scipy_operators():
    """NumPy and SciPy sparse operators infer the SciPy solver backend."""
    assert helpers._infer_backend(np.eye(2)) == "scipy"
    assert helpers._infer_backend(sp.eye(2, format="csr")) == "scipy"


def test_infer_backend_rejects_unknown_operator_types():
    """Unrecognized operator types require an explicit backend selection."""
    with pytest.raises(TypeError, match="Could not infer"):
        helpers._infer_backend(object())
