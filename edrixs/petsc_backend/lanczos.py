import numpy as np
from scipy.sparse.linalg import aslinearoperator
from petsc4py import PETSc
from slepc4py import SLEPc

def lanczos_tridiagonal_scipy(H, v0, m):
    r"""
    Perform the Lanczos tridiagonalization of a Hermitian operator :math:`H`
    starting from an initial vector :math:`v_0`, producing the diagonal
    coefficients :math:`\{\alpha_j\}` and off-diagonal coefficients
    :math:`\{\beta_j\}` of the associated Lanczos tridiagonal matrix.

    The Lanczos algorithm constructs an orthonormal Krylov basis
    :math:`\{ v_0, v_1, \ldots, v_{m-1} \}` defined by
    
    .. math::

        \mathcal{K}_m(H, v_0) =
        \mathrm{span}\{ v_0,\; H v_0,\; H^2 v_0,\; \ldots,\; H^{m-1} v_0 \}.

    In this basis, the projection of :math:`H` takes the tridiagonal form

    .. math::

        T_m =
        \begin{pmatrix}
            \alpha_0 & \beta_0  & 0        & \cdots \\
            \beta_0  & \alpha_1 & \beta_1  & \cdots \\
            0        & \beta_1  & \alpha_2 & \cdots \\
            \vdots   & \vdots   & \vdots   & \ddots
        \end{pmatrix},

    Parameters
    ----------
    H : (n,n) array-like or LinearOperator
        Hermitian operator for which the Lanczos projection is constructed.
    v0 : (n,) array-like
        Initial seed vector :math:`v_0`, which will be normalized internally.
    m : int
        Maximum number of Lanczos iterations (Krylov dimension).

    Returns
    -------
    alphas : (k,) ndarray
        Real diagonal coefficients :math:`\alpha_j` of the Lanczos tridiagonal matrix.
    betas : (k-1,) ndarray
        Real off-diagonal coefficients :math:`\beta_j`.  The length may be
        smaller than :math:`m-1` if a lucky breakdown occurs.
    """
    H = aslinearoperator(H)
    norm_psi = np.linalg.norm(v0)
    n = H.shape[0]
    v = v0 / norm_psi
    
    alphas = np.zeros(m, dtype=float)
    betas  = np.zeros(m-1, dtype=float)

    w = H @ v
    alphas[0] = np.vdot(v, w).real
    w = w - alphas[0] * v

    neff = 1
    for j in range(1, m):
        beta = np.linalg.norm(w)
        if beta == 0:
            # lucky breakdown: actual Krylov dimension < m
            return alphas[:j], betas[:j-1]
        betas[j-1] = beta
        v_old = v
        v = w / beta
        w = H @ v - beta * v_old
        alphas[j] = np.vdot(v, w).real
        w = w - alphas[j] * v
        neff += 1

    return alphas[:neff], betas[:neff-1], norm_psi**2


def lanczos_tridiagonal(H, v, nkryl):
    r"""
    Perform the Lanczos tridiagonalization of a Hermitian operator :math:`H`
    starting from an initial vector :math:`v_0`, producing the diagonal
    coefficients :math:`\{\alpha_j\}` and off-diagonal coefficients
    :math:`\{\beta_j\}` of the associated Lanczos tridiagonal matrix.

    The Lanczos algorithm constructs an orthonormal Krylov basis
    :math:`\{ v_0, v_1, \ldots, v_{m-1} \}` defined by
    
    .. math::

        \mathcal{K}_m(H, v_0) =
        \mathrm{span}\{ v_0,\; H v_0,\; H^2 v_0,\; \ldots,\; H^{m-1} v_0 \}.

    In this basis, the projection of :math:`H` takes the tridiagonal form

    .. math::

        T_m =
        \begin{pmatrix}
            \alpha_0 & \beta_0  & 0        & \cdots \\
            \beta_0  & \alpha_1 & \beta_1  & \cdots \\
            0        & \beta_1  & \alpha_2 & \cdots \\
            \vdots   & \vdots   & \vdots   & \ddots
        \end{pmatrix},

    Parameters
    ----------
    H : (n,n) PETSc sparse matrix
        Hermitian operator for which the Lanczos projection is constructed.
    v0 : (n,) PETSc sparse vector
        Initial seed vector :math:`v_0`, which will be normalized internally.
    nkryl : int
        Maximum number of Lanczos iterations (Krylov dimension).

    Returns
    -------
    alphas : (k,) ndarray
        Real diagonal coefficients :math:`\alpha_j` of the Lanczos tridiagonal matrix.
    betas : (k-1,) ndarray
        Real off-diagonal coefficients :math:`\beta_j`.  The length may be
        smaller than :math:`m-1` if a lucky breakdown occurs.
    norm**2: float
        normalization factor. 
    """
    norm = v.normalize()
    
    alphas = np.zeros(nkryl, dtype=float)
    betas  = np.zeros(nkryl - 1, dtype=float)
    
    w = H.createVecLeft()
    H.mult(v, w)
    alphas[0] = v.dot(w).real
    
    w.axpy(-alphas[0], v)
    
    neff = 1
    for j in range(1, nkryl):
        beta = w.norm()
        if beta == 0:
            # lucky breakdown: actual Krylov dimension < m
            return alphas[:j], betas[:j-1], norm**2
        betas[j-1] = beta
        v_old = v.copy()
        v = w.copy()
        v.scale(1./beta)
        w = H @ v - beta * v_old
        alphas[j] = v.dot(w).real
        w.axpy(-alphas[j], v)
        neff += 1
    return alphas, betas, norm**2
