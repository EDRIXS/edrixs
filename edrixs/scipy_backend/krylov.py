# -*- coding: utf-8 -*-
"""SciPy-based Krylov/Lanczos utilities.

This module provides low-level routines for evaluating projected resolvents
and spectral functions from a Lanczos tridiagonal representation. The
functions are independent of any specific spectroscopy solver, but use
SciPy's linear-operator and banded-solver interfaces. High-level SciPy backend
routines can build problem-specific seed vectors and then use these utilities
to evaluate the corresponding dynamical response.

The central convention used here is

    G(z) = <psi | (z I - H)^(-1) | psi>,

where ``psi`` is the original, generally unnormalized, Lanczos seed vector.
With ``z = omega + E0 + i * eta``, the corresponding broadened spectral
function is

    A(omega) = -Im G(omega + E0 + i eta) / pi.
"""

import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse.linalg import aslinearoperator

__all__ = [
    "lanczos_tridiagonal",
    "resolvent_from_tridiag",
    "spectral_function_from_tridiag",
]


def lanczos_tridiagonal(H, v0, m):
    """Construct a Lanczos tridiagonal representation of a Hermitian operator.

    Starting from a nonzero seed vector ``v0``, this routine builds an
    orthonormal Krylov basis for

    ``span{v0, H v0, H**2 v0, ...}``

    and returns the diagonal and off-diagonal elements of the tridiagonal
    projection of ``H`` in that basis.  The seed vector is normalized
    internally, but the squared norm of the original seed is returned so that
    projected spectral weights can be reconstructed.

    Parameters
    ----------
    H : array_like, sparse matrix, or scipy.sparse.linalg.LinearOperator
        Square Hermitian operator to tridiagonalize.  Only matrix-vector
        products with ``H`` are required.
    v0 : array_like, shape (n,)
        Nonzero initial seed vector.  This may be real or complex and is
        normalized internally.
    m : int
        Maximum Krylov dimension.  The actual number of Lanczos vectors may be
        smaller if a lucky breakdown occurs.  Values larger than the Hilbert
        space dimension are clipped to ``H.shape[0]``.

    Returns
    -------
    alphas : ndarray, shape (k,)
        Diagonal entries of the Lanczos tridiagonal matrix.
    betas : ndarray, shape (k - 1,)
        Off-diagonal entries of the Lanczos tridiagonal matrix.
    norm : float
        Squared norm of the original seed vector, ``||v0||**2``.

    Raises
    ------
    ValueError
        If ``H`` is not square, ``m < 1``, ``v0`` is not one-dimensional, the
        dimensions of ``H`` and ``v0`` are incompatible, or ``v0`` has zero
        norm.

    Notes
    -----
    The returned coefficients define the tridiagonal matrix

    ``T = diag(alphas) + diag(betas, 1) + diag(betas, -1)``.

    For a full Krylov space, the projected resolvent constructed from ``T`` is
    equivalent to the exact projected resolvent

    ``<v0 | (z I - H)^(-1) | v0>``.

    This implementation does not perform explicit reorthogonalization, which
    is often sufficient for moderate Krylov dimensions but may need to be
    revisited for very large calculations or nearly degenerate spectra.
    """
    H = aslinearoperator(H)

    if H.shape[0] != H.shape[1]:
        raise ValueError("H must be a square operator")

    m = int(m)
    if m < 1:
        raise ValueError("m must be a positive integer")

    v0 = np.asarray(v0)
    if v0.ndim != 1:
        raise ValueError("v0 must be a one-dimensional array")
    if v0.shape[0] != H.shape[1]:
        raise ValueError(
            "dimension mismatch: H has shape {} but v0 has shape {}".format(
                H.shape, v0.shape
            )
        )

    norm = np.linalg.norm(v0)
    if norm == 0:
        raise ValueError("v0 must have nonzero norm")

    m = min(m, H.shape[0])

    v = v0.astype(np.result_type(v0.dtype, np.complex128), copy=False) / norm
    v_old = np.zeros_like(v)

    alphas = np.zeros(m, dtype=float)
    betas = np.zeros(max(m - 1, 0), dtype=float)

    for j in range(m):
        w = H @ v
        if j > 0:
            w = w - betas[j - 1] * v_old

        alpha = np.vdot(v, w).real
        alphas[j] = alpha
        w = w - alpha * v

        if j == m - 1:
            return alphas, betas, norm**2

        beta = np.linalg.norm(w)
        if beta == 0:
            return alphas[: j + 1], betas[:j], norm**2

        betas[j] = beta
        v_old = v
        v = w / beta

    return alphas, betas, norm**2


def _validate_tridiag(alphas, betas):
    """Validate and normalize Lanczos tridiagonal coefficient arrays.

    Parameters
    ----------
    alphas : array_like, shape (m,)
        Diagonal entries of the tridiagonal matrix.
    betas : array_like, shape (m - 1,)
        Off-diagonal entries of the tridiagonal matrix.

    Returns
    -------
    alphas : ndarray, shape (m,)
        Validated one-dimensional float array.
    betas : ndarray, shape (m - 1,)
        Validated one-dimensional float array.

    Raises
    ------
    ValueError
        If either input is not one-dimensional, ``alphas`` is empty, or
        ``len(betas) != len(alphas) - 1``.
    """
    alphas = np.asarray(alphas, dtype=float)
    betas = np.asarray(betas, dtype=float)

    if alphas.ndim != 1:
        raise ValueError("alphas must be one-dimensional")
    if betas.ndim != 1:
        raise ValueError("betas must be one-dimensional")
    if len(alphas) == 0:
        raise ValueError("alphas must contain at least one element")
    if len(betas) != max(len(alphas) - 1, 0):
        raise ValueError("betas must have length len(alphas) - 1")

    return alphas, betas


def resolvent_from_tridiag(alphas, betas, z, norm=1.0):
    """Evaluate a projected resolvent from Lanczos tridiagonal coefficients.

    This routine computes

    ``G(z) = norm * [(z I - T)^(-1)]_00``,

    where ``T`` is the tridiagonal Lanczos matrix defined by ``alphas`` and
    ``betas``.  If ``alphas`` and ``betas`` were generated by
    :func:`lanczos_tridiagonal` using a seed vector ``psi`` and ``norm`` is the
    returned squared norm ``||psi||**2``, then ``G(z)`` approximates

    ``<psi | (z I - H)^(-1) | psi>``.

    Parameters
    ----------
    alphas : array_like, shape (m,)
        Diagonal entries of the Lanczos tridiagonal matrix.
    betas : array_like, shape (m - 1,)
        Off-diagonal entries of the Lanczos tridiagonal matrix.
    z : complex or array_like of complex
        Complex energy argument or array of energy arguments.  The convention
        is ``z I - H``.  For retarded response functions one typically uses
        ``z = omega + E0 + 1j * eta`` with ``eta > 0``.
    norm : float or complex, optional
        Overall spectral weight.  For coefficients generated by
        :func:`lanczos_tridiagonal`, this should usually be the returned
        ``norm`` value, equal to the squared norm of the original seed vector.
        The default is ``1.0``.

    Returns
    -------
    G : complex or ndarray of complex
        Projected resolvent evaluated at ``z``.  A scalar is returned for
        scalar input; otherwise the output has the same shape as ``z``.

    Raises
    ------
    ValueError
        If the tridiagonal coefficient arrays have incompatible shapes.

    Notes
    -----
    Each value of ``z`` requires solving only a small tridiagonal linear system
    of size ``len(alphas)``.  The full many-body Hamiltonian is not used by
    this function.
    """
    alphas, betas = _validate_tridiag(alphas, betas)

    z_arr = np.asarray(z, dtype=complex)
    scalar_input = z_arr.ndim == 0
    z_flat = z_arr.reshape(-1)

    m = len(alphas)
    e1 = np.zeros(m, dtype=complex)
    e1[0] = 1.0

    out = np.empty(z_flat.shape, dtype=complex)

    for iz, zi in enumerate(z_flat):
        ab = np.zeros((3, m), dtype=complex)
        ab[1, :] = zi - alphas
        ab[0, 1:] = -betas
        ab[2, :-1] = -betas

        y = solve_banded((1, 1), ab, e1)
        out[iz] = norm * y[0]

    if scalar_input:
        return out[0]
    return out.reshape(z_arr.shape)


def spectral_function_from_tridiag(alphas, betas, omega, eta, e0=0.0, norm=1.0):
    """Evaluate a broadened spectral function from Lanczos coefficients.

    This routine evaluates the retarded spectral function

    ``A(omega) = -Im G(omega + e0 + i eta) / pi``,

    where ``G`` is computed by :func:`resolvent_from_tridiag`.  With
    coefficients generated from a seed vector ``psi = O |g>`` and
    ``e0 = E_g``, this corresponds to the dynamical correlation function

    ``A(omega) = -Im <psi | (omega + E_g + i eta - H)^(-1) | psi> / pi``.

    Parameters
    ----------
    alphas : array_like, shape (m,)
        Diagonal entries of the Lanczos tridiagonal matrix.
    betas : array_like, shape (m - 1,)
        Off-diagonal entries of the Lanczos tridiagonal matrix.
    omega : float or array_like of float
        Real frequency or frequency grid at which to evaluate the spectral
        function.
    eta : float
        Positive Lorentzian broadening parameter.
    e0 : float, optional
        Energy shift, typically the energy of the initial state.  The default
        is ``0.0``.
    norm : float or complex, optional
        Overall spectral weight.  For coefficients generated by
        :func:`lanczos_tridiagonal`, this should usually be the returned
        squared seed norm.  The default is ``1.0``.

    Returns
    -------
    spectrum : float or ndarray of float
        Broadened spectral function evaluated on ``omega``.  A scalar is
        returned for scalar input; otherwise the output has the same shape as
        ``omega``.

    Raises
    ------
    ValueError
        If ``eta`` is not positive or the tridiagonal coefficient arrays have
        incompatible shapes.

    Notes
    -----
    The sign convention is consistent with a retarded resolvent
    ``(omega + e0 + i eta - H)^(-1)``.  The output is positive for positive
    spectral weights when ``H`` is Hermitian and ``eta > 0``.
    """
    if eta <= 0:
        raise ValueError("eta must be positive")

    z = np.asarray(omega, dtype=float) + e0 + 1j * eta
    g = resolvent_from_tridiag(alphas, betas, z, norm=norm)
    return -np.imag(g) / np.pi
