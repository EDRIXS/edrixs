"""Unit tests for the response-building stages inside SciPy XAS and RIXS.

After photon transitions create a start vector, the SciPy workflows compress
Hamiltonian action into a short Krylov representation.  These tests isolate
that numerical stage before broadening and final spectrum assembly occur.
"""

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
    """Check the intermediate-response stage used by ``xas_krylov_scipy``.

    XAS applies a photon transition to an initial state and then uses a Krylov
    representation of the intermediate Hamiltonian to evaluate its response
    across incident energies.  This test checks that response against a direct
    dense solve before it is converted into an XAS spectrum.
    """
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
    """Check early completion of the Krylov stage used by XAS and RIXS.

    When the transition-generated vector is already an eigenvector, no further
    Hamiltonian directions are needed.  This test ensures the response stage
    stops cleanly rather than creating spurious poles in the final spectrum.
    """
    hmat = np.diag([1.0, 2.0, 4.0])
    seed = np.array([0.0, 3.0, 0.0])

    alpha, beta, norm = lanczos_tridiagonal(hmat, seed, m=3)

    assert_allclose(alpha, [2.0])
    assert beta.size == 0
    assert norm == pytest.approx(9.0)


def test_spectral_function_has_scalar_and_array_contracts():
    """Check the pole-evaluation utility used during spectrum assembly.

    The Krylov stages return compact pole information, and the public spectral
    workflows evaluate it either at one energy or over an energy grid.  This
    test checks both input forms before the values are placed in XAS/RIXS arrays.
    """
    alpha = np.array([1.5])
    beta = np.array([])

    scalar = spectral_function_from_tridiag(alpha, beta, 1.5, eta=0.2)
    array = spectral_function_from_tridiag(alpha, beta, [1.4, 1.5], eta=0.2)

    assert np.isscalar(scalar)
    assert array.shape == (2,)
    assert scalar > 0


@pytest.mark.parametrize(
    "hmat, seed, msteps, message",
    [
        (np.ones((2, 3)), np.ones(3), 2, "square"),
        (np.eye(2), np.ones(2), 0, "positive"),
        (np.eye(2), np.ones((2, 1)), 2, "one-dimensional"),
        (np.eye(2), np.ones(3), 2, "dimension mismatch"),
        (np.eye(2), np.zeros(2), 2, "nonzero norm"),
    ],
)
def test_lanczos_validation(hmat, seed, msteps, message):
    """Reject invalid inputs before the XAS/RIXS Krylov stage begins.

    These checks protect the internal response calculation from malformed
    Hamiltonians, start vectors, and iteration limits that would otherwise
    produce obscure failures while constructing the final spectrum.
    """
    with pytest.raises(ValueError, match=message):
        lanczos_tridiagonal(hmat, seed, msteps)


@pytest.mark.parametrize(
    "alpha, beta",
    [
        ([], []),
        ([[1.0]], []),
        ([1.0], [[0.1]]),
        ([1.0, 2.0], []),
    ],
)
def test_resolvent_rejects_invalid_tridiagonal_shapes(alpha, beta):
    """Reject malformed pole data before an XAS response is evaluated.

    The response evaluator receives diagonal and off-diagonal coefficients
    from the Krylov stage.  This test ensures inconsistent shapes are stopped
    before they can be interpreted as physical spectral weight.
    """
    with pytest.raises(ValueError):
        resolvent_from_tridiag(alpha, beta, 1j)


def test_spectral_function_requires_positive_eta():
    """Require physical broadening when evaluating Krylov spectral poles.

    Pole data are broadened to form a finite spectrum.  This test ensures the
    internal evaluator rejects a nonpositive width rather than returning an
    undefined value to the XAS or RIXS spectrum assembly stage.
    """
    with pytest.raises(ValueError, match="positive"):
        spectral_function_from_tridiag([1.0], [], 0.0, eta=0.0)
