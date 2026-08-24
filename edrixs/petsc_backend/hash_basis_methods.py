"""Compatibility names for the unimplemented PETSc backend.

Fock-basis representations and ranking are backend independent and live in
:mod:`edrixs.fock_basis`.  ``scipy_edrixs`` intentionally keeps PETSc matrix
construction as a stub, so this module exposes only the historical basis names
without introducing a second PETSc implementation.
"""

from __future__ import annotations

from ..fock_basis import FockBinByN, get_fock_basis_combinadic

__all__ = [
    'FockBinByN', 'get_fock_basis_petsc', 'build_op_petsc_matrix',
]


# Backwards-compatible name used by the PETSc branch.  The representation
# itself is backend independent and is shared with SciPy.
get_fock_basis_petsc = get_fock_basis_combinadic


def build_op_petsc_matrix(*args, **kwargs):
    """Raise the standard error for unavailable PETSc matrix construction."""
    raise NotImplementedError(
        "The PETSc backend contract is present, but build_op_petsc has not yet "
        "been implemented"
    )
