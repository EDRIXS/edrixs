"""Unit tests for backend-neutral helpers used by :mod:`edrixs.solvers`."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

import edrixs._solvers_helpers as helpers


@pytest.mark.parametrize("scale", [1.0, 1e-12])
def test_dense_to_sparse_u_preserves_flattening_convention(scale):
    """Dense and flattened sparse Coulomb data must use the same ordering."""
    rng = np.random.default_rng(3)
    umat = rng.normal(size=(3, 3, 3, 3))
    umat[np.abs(umat) < 0.6] = 0
    umat = umat.astype(complex) * scale * (1 + 1j)

    actual = helpers._umat_dense_to_sparse(umat).toarray()

    np.testing.assert_array_equal(actual, umat.reshape(9, 9))


@pytest.mark.parametrize("scale", [1.0, 1e-12])
def test_sparse_siam_embedding_matches_dense_embedding(scale):
    """Sparse SIAM embedding must match the dense reference placement."""
    rng = np.random.default_rng(8)
    compact = rng.normal(size=(4, 4, 4, 4))
    compact[np.abs(compact) < 0.6] = 0
    compact = compact.astype(complex) * scale * (1 + 1j)
    dense = helpers._embed_impurity_core_umat(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )
    sparse = helpers._embed_impurity_core_umat_sparse(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )

    np.testing.assert_array_equal(sparse.toarray(), dense.reshape(64, 64))


def test_expand_broadening_accepts_scalar_or_exact_length_array():
    """Broadening normalization accepts a scalar or one value per grid point."""
    assert_allclose(helpers._expand_broadening(0.2, 3, "gamma"), [0.2] * 3)
    assert_allclose(
        helpers._expand_broadening([0.1, 0.2], 2, "gamma"),
        [0.1, 0.2],
    )

    with pytest.raises(ValueError, match="shape"):
        helpers._expand_broadening([0.1], 2, "gamma")


def test_infer_backend_recognizes_dense_and_scipy_operators():
    """NumPy and SciPy sparse operators infer their respective backends."""
    assert helpers._infer_backend(np.eye(2)) == "dense"
    assert helpers._infer_backend(sp.eye(2, format="csr")) == "scipy"


def test_infer_backend_rejects_unknown_operator_types():
    """Unrecognized operator types require an explicit backend selection."""
    with pytest.raises(TypeError, match="Could not infer"):
        helpers._infer_backend(object())
