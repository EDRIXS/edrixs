"""SciPy-specific many-body operator construction helpers.

Fock-basis representations and ranking live in :mod:`edrixs.fock_basis`.
This module owns SciPy sparse-matrix materialization.  The shared operator code
is used only as a matrix-entry kernel over backend-selected basis ranges.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .._operator_builder import prepare_operator_entry_kernel

__all__ = ['build_op_scipy_matrix']


def _assemble_csr(rows, cols, data, shape):
    """Materialize the final SciPy sparse operator using the current COO path."""
    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=shape,
        dtype=np.complex128,
    ).tocsr()


def build_op_scipy_matrix(
        emat, umat, lb, rb=None, *, tol=1e-10, use_numba=False):
    """Build the final SciPy CSR matrix using the SciPy assembly policy."""
    if rb is None:
        rb = lb

    build_range = prepare_operator_entry_kernel(
        emat,
        umat,
        lb,
        rb,
        tol_e=tol,
        tol_u=tol,
        use_numba=use_numba,
    )

    # SciPy currently materializes the complete COO contribution arrays before
    # converting them to CSR.  This policy is intentionally backend-local.
    rows, cols, data = build_range(0, len(rb))
    return _assemble_csr(rows, cols, data, (len(lb), len(rb)))
