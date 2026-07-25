"""Unit tests for the initial-state solution stage of the SciPy workflow.

``ed_krylov_scipy`` produces the low-energy initial states subsequently passed
to ``xas_krylov_scipy`` and ``rixs_krylov_scipy``.  These tests isolate that
stage from model setup, transition construction, and final spectra.
"""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import ed_krylov_scipy


def test_ed_krylov_returns_lowest_sorted_eigenpairs():
    """Check the initial-state eigenpairs supplied to XAS and RIXS.

    The spectral workflows require the lowest retained energies and matching
    eigenvectors in ascending order.  This test verifies that contract before
    those states are weighted, acted on by photons, or used in spectra.
    """
    diagonal = np.array([-3.0, -1.0, 0.5, 1.0, 2.0, 4.0, 5.0, 8.0])
    hmat = sp.diags(diagonal, format="csr")

    evals, evecs = ed_krylov_scipy(
        hmat,
        num_gs=2,
        blocksize=2,
        initial_guess=np.eye(len(diagonal), 2),
        tol=1e-12,
    )

    assert_allclose(evals, diagonal[:2], atol=1e-12)
    assert_allclose(hmat @ evecs, evecs * evals, atol=1e-11)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"num_gs": 0}, "positive"),
        ({"num_gs": 2, "blocksize": 1}, "greater than or equal"),
        ({"num_gs": 1, "blocksize": 9}, "cannot exceed"),
        (
            {"num_gs": 1, "blocksize": 1, "initial_guess": np.ones((7, 1))},
            "initial_guess",
        ),
    ],
)
def test_ed_krylov_validation(kwargs, message):
    """Reject invalid ED settings before initial states are calculated.

    These public-boundary checks stop impossible state counts, block sizes, and
    initial guesses before the eigensolver runs and before invalid states can
    propagate into XAS or RIXS calculations.
    """
    with pytest.raises(ValueError, match=message):
        ed_krylov_scipy(np.eye(8), **kwargs)
