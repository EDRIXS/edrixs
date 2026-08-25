"""PETSc-specific many-body operator construction helpers.

Fock-basis representations and ranking live in :mod:`edrixs.fock_basis`.
This module owns PETSc work decomposition, preallocation, insertion, and final
matrix assembly.  The shared operator code is used only as a bounded-range
matrix-entry kernel.
"""

from __future__ import annotations

import numpy as np

from .._operator_builder import prepare_operator_entry_kernel
from ..fock_basis import FockBinByN, get_fock_basis_combinadic

__all__ = [
    'FockBinByN', 'get_fock_basis_petsc', 'build_op_petsc_matrix',
]


# Backwards-compatible name from the initial PETSc implementation.
get_fock_basis_petsc = get_fock_basis_combinadic


def _default_preallocation(emat, umat, nr, tol_e, tol_u):
    """Return today's simple PETSc preallocation policy.

    Keeping this policy separate makes it possible to replace the heuristic
    with exact local diagonal/off-diagonal counts without touching the shared
    operator kernels.
    """
    ne = int(np.count_nonzero(np.abs(emat) > tol_e)) if emat is not None else 0
    if umat is None:
        nu = 0
    elif hasattr(umat, 'data') and hasattr(umat, 'tocoo'):
        nu = int(np.count_nonzero(np.abs(umat.data) > tol_u))
    else:
        nu = int(np.count_nonzero(np.abs(umat) > tol_u))
    return max(8, min(nr, 16 + ne + 2 * min(nu, 32)))


def _create_matrix(PETSc, comm, nl, nr, preallocation):
    """Create the final PETSc matrix using the current CPU-AIJ policy."""
    hmat = PETSc.Mat().create(comm=comm)
    hmat.setSizes(((None, nl), (None, nr)))
    hmat.setType(PETSc.Mat.Type.AIJ)
    hmat.setPreallocationNNZ(preallocation)
    hmat.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    hmat.setUp()
    return hmat


def _work_ranges(hmat, chunk_cols):
    """Yield today's MPI work decomposition as bounded owned-column ranges."""
    cstart, cend = hmat.getOwnershipRangeColumn()
    for start in range(cstart, cend, chunk_cols):
        yield start, min(start + chunk_cols, cend)


def _insert_entries(hmat, rows, cols, vals, PETSc):
    """Insert one bounded contribution buffer using PETSc ``setValues``."""
    if rows.size == 0:
        return

    rows = np.asarray(rows, dtype=PETSc.IntType)
    cols = np.asarray(cols, dtype=PETSc.IntType)
    vals = np.asarray(vals, dtype=PETSc.ScalarType)

    # AIJ is row-oriented.  Sorting is bounded by one work chunk rather than
    # the complete local sparse matrix.  A future GPU/COO insertion strategy
    # can replace this helper without changing the operator kernels.
    order = np.argsort(rows, kind='mergesort')
    rows = rows[order]
    cols = cols[order]
    vals = vals[order]

    start = 0
    while start < rows.size:
        row = int(rows[start])
        end = start + 1
        while end < rows.size and rows[end] == rows[start]:
            end += 1
        hmat.setValues(
            row,
            cols[start:end],
            vals[start:end],
            addv=PETSc.InsertMode.ADD_VALUES,
        )
        start = end


def _finalize_matrix(hmat, mat_type):
    """Finalize PETSc assembly using the current post-assembly conversion."""
    hmat.assemblyBegin()
    hmat.assemblyEnd()
    if mat_type is not None:
        hmat = hmat.convert(mat_type)
    return hmat


def build_op_petsc_matrix(
    emat,
    umat,
    lb,
    rb=None,
    *,
    comm,
    tol_e=1e-10,
    tol_u=1e-10,
    nnz_guess_per_row=None,
    mat_type=None,
    assembly_chunk_cols=4096,
    use_numba=False,
):
    """Build the final PETSc matrix from backend-selected bounded work units."""
    from petsc4py import PETSc

    if rb is None:
        rb = lb
    if lb.norbs != rb.norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    assembly_chunk_cols = int(assembly_chunk_cols)
    if assembly_chunk_cols < 1:
        raise ValueError("assembly_chunk_cols must be a positive integer")

    nl, nr = len(lb), len(rb)
    if nnz_guess_per_row is None:
        nnz_guess_per_row = _default_preallocation(
            emat, umat, nr, tol_e, tol_u
        )

    # PETSc owns the matrix before any many-body matrix entries are generated.
    hmat = _create_matrix(
        PETSc, comm, nl, nr, nnz_guess_per_row
    )

    build_range = prepare_operator_entry_kernel(
        emat,
        umat,
        lb,
        rb,
        tol_e=tol_e,
        tol_u=tol_u,
        use_numba=use_numba,
    )

    for start, end in _work_ranges(hmat, assembly_chunk_cols):
        rows, cols, vals = build_range(start, end)
        _insert_entries(hmat, rows, cols, vals, PETSc)

    return _finalize_matrix(hmat, mat_type)
