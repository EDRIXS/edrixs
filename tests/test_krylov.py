"""Unit tests for the backend-independent Krylov response utilities."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse import csr_matrix

from edrixs.krylov import (
    lanczos_tridiagonal,
    resolvent_from_tridiag,
    spectral_function_from_tridiag,
)


def test_full_lanczos_resolvent_matches_direct_solve():
    """A full Krylov representation reproduces a direct resolvent solve."""
    rng = np.random.default_rng(12)
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    hmat = (raw + raw.conj().T) / 2
    seed = rng.normal(size=4) + 1j * rng.normal(size=4)
    zvals = np.array([-0.4 + 0.3j, 0.2 + 0.15j, 1.1 + 0.5j])

    alpha, beta, norm = lanczos_tridiagonal(csr_matrix(hmat), seed, m=4)
    actual = resolvent_from_tridiag(alpha, beta, zvals, norm=norm)
    expected = np.array(
        [
            np.vdot(seed, np.linalg.solve(zval * np.eye(4) - hmat, seed))
            for zval in zvals
        ]
    )

    assert_allclose(actual, expected, rtol=2e-11, atol=2e-11)


def test_lanczos_lucky_breakdown_for_eigenvector_seed():
    """Lanczos stops cleanly when its seed is already an eigenvector."""
    hmat = np.diag([1.0, 2.0, 4.0])
    seed = np.array([0.0, 3.0, 0.0])

    alpha, beta, norm = lanczos_tridiagonal(hmat, seed, m=3)

    assert_allclose(alpha, [2.0])
    assert beta.size == 0
    assert norm == pytest.approx(9.0)


def test_spectral_function_has_scalar_and_array_contracts():
    """Pole evaluation preserves scalar and energy-grid return contracts."""
    alpha = np.array([1.5])
    beta = np.array([])

    scalar = spectral_function_from_tridiag(alpha, beta, 1.5, eta=0.2)
    array = spectral_function_from_tridiag(
        alpha, beta, [1.4, 1.5], eta=0.2
    )

    assert np.isscalar(scalar)
    assert array.shape == (2,)
    assert scalar > 0


@pytest.mark.parametrize(
    ("hmat", "seed", "msteps", "message"),
    [
        (np.ones((2, 3)), np.ones(3), 2, "square"),
        (np.eye(2), np.ones(2), 0, "positive"),
        (np.eye(2), np.ones((2, 1)), 2, "one-dimensional"),
        (np.eye(2), np.ones(3), 2, "dimension mismatch"),
        (np.eye(2), np.zeros(2), 2, "nonzero norm"),
    ],
)
def test_lanczos_validation(hmat, seed, msteps, message):
    """Malformed Krylov inputs are rejected with useful errors."""
    with pytest.raises(ValueError, match=message):
        lanczos_tridiagonal(hmat, seed, msteps)


@pytest.mark.parametrize(
    ("alpha", "beta"),
    [
        ([], []),
        ([[1.0]], []),
        ([1.0], [[0.1]]),
        ([1.0, 2.0], []),
    ],
)
def test_resolvent_rejects_invalid_tridiagonal_shapes(alpha, beta):
    """Malformed pole arrays are rejected before response evaluation."""
    with pytest.raises(ValueError):
        resolvent_from_tridiag(alpha, beta, 1j)


def test_spectral_function_requires_positive_eta():
    """Spectrum evaluation requires a strictly positive broadening."""
    with pytest.raises(ValueError, match="positive"):
        spectral_function_from_tridiag([1.0], [], 0.0, eta=0.0)
