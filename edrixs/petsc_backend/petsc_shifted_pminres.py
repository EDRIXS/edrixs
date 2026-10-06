"""Reusable PETSc solver for a Hermitian operator plus a scalar imaginary shift."""

import cmath
import ctypes as ct

import numpy as np
from petsc4py import PETSc


class ShiftedPMINRES:
    """PETSc Python KSP context for ``(H + i*gamma*I) x = b``.

    ``H`` must be Hermitian and ``gamma`` real. The caller supplies the
    already shifted matrix through ``ksp.setOperators``. Only PCNONE and a
    zero initial guess are supported; Hermiticity is a caller responsibility.
    Attach one instance to one KSP with ``setPythonContext`` and reuse that
    KSP for successive matrices and right-hand sides of the same layout.
    The Python recurrence retains five work vectors until reset or destroy.

    Parameters
    ----------
    library : str or pathlib.Path, optional
        Shared library exporting ``shifted_pminres``. If omitted, use the
        Python recurrence. The native ABI is ``int(Mat, Vec, Vec, double,
        PetscInt, PetscInt*, PetscInt*)`` (operator, RHS, solution, absolute
        tolerance, iteration limit, iterations, convergence reason). It must
        link to the same PETSc build as petsc4py, with double-precision PetscReal
        and matching scalar type and PetscInt width. Native work allocation is library-owned.
        The native implementation requires ``rtol=0``.

    Notes
    -----
    Solves and cleanup are collective on the KSP communicator. The supplied
    matrix and RHS are borrowed and never destroyed or modified. The solution
    is overwritten. Convergence is checked using the true residual norm;
    iteration count, residual norm and reason are available on the KSP.
    Custom KSP monitors and convergence callbacks are not invoked.
    """

    def __init__(self, library=None):
        """Initialize an empty workspace and optionally load the native kernel."""
        self._vectors = []
        self._layout = None
        self._function = None
        if library is not None:
            if np.dtype(PETSc.RealType) != np.dtype('float64'):
                raise ValueError('The native PMINRES ABI requires double precision')
            self._integer = ct.c_int32 if np.dtype(PETSc.IntType).itemsize == 4 else ct.c_int64
            self._library = ct.CDLL(str(library))
            self._function = self._library.shifted_pminres
            self._function.argtypes = [ct.c_void_p] * 3 + [
                ct.c_double, self._integer,
                ct.POINTER(self._integer), ct.POINTER(self._integer),
            ]
            self._function.restype = ct.c_int

    def reset(self, ksp=None):
        """Release cached work vectors; the next solve allocates a new workspace."""
        for vector in self._vectors:
            vector.destroy()
        self._vectors = []
        self._layout = None

    def destroy(self, ksp=None):
        """Release owned PETSc resources when the KSP context is destroyed."""
        self.reset(ksp)

    def _workspace(self, b):
        """Return work vectors matching the RHS layout and vector implementation."""
        layout = (b.getSize(), b.getLocalSize(), b.getOwnershipRange(), b.getType())
        if layout != self._layout:
            self.reset()
            for _ in range(1 if self._function is not None else 5):
                self._vectors.append(b.duplicate())
            self._layout = layout
        return self._vectors

    def solve(self, ksp, b, x):
        """Solve the current KSP system into ``x`` and publish convergence status.

        Parameters
        ----------
        ksp : petsc4py.PETSc.KSP
            Owning Python KSP supplying the operator and tolerances.
        b, x : petsc4py.PETSc.Vec
            Right-hand side and distinct, layout-compatible output vector.

        The Python path accepts absolute and relative tolerances, using
        ``max(atol, rtol * ||b||)``. A failed solve sets a negative KSP reason.
        Invalid preconditioners or initial guesses raise ``ValueError``.
        """
        if ksp.getPC().getType() != 'none' or ksp.getInitialGuessNonzero():
            raise ValueError('ShiftedPMINRES needs PCNONE and a zero initial guess')
        matrix, _ = ksp.getOperators()
        rtol, atol, _, maxits = ksp.getTolerances()
        if self._function is not None:
            if rtol != 0:
                raise ValueError('Native ShiftedPMINRES requires an absolute tolerance')
            residual, = self._workspace(b)
            iterations, reason = self._integer(), self._integer()
            x.set(0)
            error = self._function(matrix.handle, b.handle, x.handle, atol, maxits,
                                   ct.byref(iterations), ct.byref(reason))
            if error:
                raise PETSc.Error(error)
            matrix.mult(x, residual)
            residual.axpy(-1, b)
            norm = residual.norm()
            status = reason.value
            if status > 0 and not norm <= atol:
                status = PETSc.KSP.ConvergedReason.DIVERGED_BREAKDOWN
            ksp.setResidualNorm(norm)
            ksp.setIterationNumber(iterations.value)
            ksp.setConvergedReason(status)
            return
        vectors = self._workspace(b)
        v, vprev, w, wprev, work = vectors
        for vec in vectors:
            vec.set(0)
        x.set(0)
        b.copy(v)
        beta = b.norm()
        threshold = max(atol, rtol * beta)
        eta = estimate = beta
        gamma0 = gamma1 = 1.0 + 0j
        sigma0 = sigma1 = 0j
        reason = PETSc.KSP.ConvergedReason.DIVERGED_MAX_IT
        iterations = 0
        if beta <= threshold:
            reason = (PETSc.KSP.ConvergedReason.CONVERGED_ATOL
                      if beta <= atol else PETSc.KSP.ConvergedReason.CONVERGED_RTOL)
        else:
            for iterations in range(1, maxits + 1):
                v.scale(1 / beta)
                matrix.mult(v, work)
                # PETSc conjugates its SECOND dot-product argument.
                alpha = work.dot(v)
                vprev.axpby(1, -beta, work)
                vprev.axpy(-alpha, v)
                beta_next = vprev.norm()
                delta = gamma0 * alpha - gamma1 * sigma0 * beta
                rho1 = cmath.sqrt(delta * delta + beta_next * beta_next)
                if abs(rho1) == 0:
                    reason = PETSc.KSP.ConvergedReason.DIVERGED_BREAKDOWN
                    break
                rho2 = sigma0 * alpha + gamma1 * gamma0 * beta
                rho3 = sigma1 * beta
                gamma1, sigma1 = gamma0, sigma0
                gamma0, sigma0 = delta / rho1, beta_next / rho1
                # Fuse scaling with addition and fold in the normalization.
                wprev.axpby(1 / rho1, -rho3 / rho1, v)
                wprev.axpy(-rho2 / rho1, w)
                x.axpy(gamma0 * eta, wprev)
                estimate *= abs(sigma0)
                eta *= -sigma0
                if estimate <= threshold or beta_next == 0:
                    matrix.mult(x, work)
                    work.axpy(-1, b)
                    norm = work.norm()
                    if norm <= threshold:
                        reason = (PETSc.KSP.ConvergedReason.CONVERGED_ATOL
                                  if norm <= atol else PETSc.KSP.ConvergedReason.CONVERGED_RTOL)
                        break
                if beta_next == 0:
                    reason = PETSc.KSP.ConvergedReason.DIVERGED_BREAKDOWN
                    break
                beta = beta_next
                v, vprev = vprev, v
                w, wprev = wprev, w
        matrix.mult(x, work)
        work.axpy(-1, b)
        ksp.setIterationNumber(iterations)
        ksp.setResidualNorm(work.norm())
        ksp.setConvergedReason(reason)
