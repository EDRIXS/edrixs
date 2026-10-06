import numpy as np


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
    v : (n,) PETSc sparse vector
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
    nkryl = min(nkryl, H.getSize()[0])
    v = v.copy()
    w = H.createVecLeft()
    v_old = v.duplicate()
    workspace = (v, w, v_old)
    try:
        norm = v.normalize()
        alphas = np.zeros(nkryl, dtype=float)
        betas = np.zeros(nkryl - 1, dtype=float)

        H.mult(v, w)
        alphas[0] = v.dot(w).real
        w.axpy(-alphas[0], v)

        for j in range(1, nkryl):
            beta = w.norm()
            if beta == 0:
                # Lucky breakdown: actual Krylov dimension < m.
                return alphas[:j], betas[:j-1], norm**2
            betas[j-1] = beta
            # The residual becomes the next basis vector; recycle the oldest
            # vector as the matrix-product destination without allocating.
            v_old, v, w = v, w, v_old
            v.scale(1. / beta)
            H.mult(v, w)
            w.axpy(-beta, v_old)
            alphas[j] = v.dot(w).real
            w.axpy(-alphas[j], v)
        return alphas, betas, norm**2
    finally:
        for vector in workspace:
            vector.destroy()
