"""PETSc backend interface stub.

This module defines the backend contract without claiming an implementation.
Importing EDRIXS does not require petsc4py; it is imported lazily only for type
recognition or when a PETSc implementation is requested.
"""

from __future__ import annotations

__all__ = [
    'owns_operator_petsc',
    'get_H_petsc', 'get_transition_operators_petsc',
    'ed_petsc', 'xas_petsc', 'rixs_petsc',
]


def _petsc_module():
    """
    Import petsc4py lazily and return its PETSc module.
    """
    try:
        from petsc4py import PETSc
    except ImportError as exc:
        raise ImportError(
            "The PETSc backend requires petsc4py and a working PETSc installation"
        ) from exc
    return PETSc


def owns_operator_petsc(operator):
    """
    Return whether ``operator`` is a petsc4py matrix or linear operator.
    """
    try:
        PETSc = _petsc_module()
    except ImportError:
        return False
    return isinstance(operator, (PETSc.Mat, PETSc.Vec))


def _not_implemented(operation):
    """
    Raise the standard PETSc-backend stub error.
    """
    raise NotImplementedError(
        "The PETSc backend contract is present, but {} has not yet been "
        "implemented".format(operation)
    )


def get_H_petsc(emat, umat, basis, *, backend_kws=None):
    """Build a PETSc many-body Hamiltonian (stub)."""
    _petsc_module()
    _not_implemented('get_H_petsc')


def get_transition_operators_petsc(
        trans_mat, basis_n, basis_i, *, backend_kws=None):
    """Build PETSc transition operators (stub)."""
    _petsc_module()
    _not_implemented('get_transition_operators_petsc')


def ed_petsc(hmat_i, num_evals=1, *, backend_kws=None):
    """Obtain low-energy eigenpairs with SLEPc/PETSc (stub)."""
    _petsc_module()
    _not_implemented('ed_petsc')


def xas_petsc(eval_i, evec_i, hmat_n, trans_op, ominc, *,
              gamma_c=0.1, thin=1.0, phi=0.0, pol_type=None,
              temperature=1.0, scatter_axis=None, backend_kws=None):
    """Calculate XAS with PETSc (stub)."""
    _petsc_module()
    _not_implemented('xas_petsc')


def rixs_petsc(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
               gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
               pol_type=None, temperature=1.0, scatter_axis=None,
               skip_gs=False, return_poles=False, backend_kws=None):
    """Calculate RIXS with PETSc (stub)."""
    _petsc_module()
    _not_implemented('rixs_petsc')
