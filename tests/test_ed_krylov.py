"""Unit tests for the public initial-state solver interface."""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import ed


def test_ed_returns_lowest_sorted_eigenpairs_and_infers_scipy():
    """The public ED interface returns matching low-energy eigenpairs."""
    diagonal = np.array([
        -3.0, -1.0, 0.5, 1.0, 2.0, 4.0,
        5.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    ])
    hmat = sp.diags(diagonal, format="csr")

    evals, evecs = ed(
        hmat,
        num_evals=2,
        backend_kws={
            "blocksize": 2,
            "initial_guess": np.eye(len(diagonal), 2),
            "tol": 1e-12,
        },
    )

    assert_allclose(evals, diagonal[:2], atol=1e-12)
    assert_allclose(hmat @ evecs, evecs * evals, atol=1e-11)


@pytest.mark.parametrize(
    ("num_evals", "backend_kws", "message"),
    [
        (0, {}, "positive"),
        (2, {"blocksize": 1}, "greater than or equal"),
        (1, {"blocksize": 9}, "cannot exceed"),
        (
            1,
            {"blocksize": 1, "initial_guess": np.ones((7, 1))},
            "initial_guess",
        ),
    ],
)
def test_ed_validation_uses_public_backend_options(
        num_evals, backend_kws, message):
    """Invalid backend-specific ED settings fail through the public interface."""
    with pytest.raises(ValueError, match=message):
        ed(
            np.eye(8),
            num_evals=num_evals,
            backend="scipy",
            backend_kws=backend_kws,
        )
