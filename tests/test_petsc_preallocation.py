"""Allocation and value regressions; also run with mpirun -n 8 pytest."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from edrixs.fock_basis import FockBasisSpec
from edrixs.solvers import build_op
from edrixs.petsc_backend.hash_basis_methods import (
    _create_matrix, _insert_entries, _finalize_matrix,
)


PETSc = pytest.importorskip('petsc4py.PETSc')
pytest.importorskip('mpi4py.MPI')


def assert_matrix_values_and_allocations(matrix, expected):
    """Compare the assembled global matrix and check for sparse storage growth."""
    lo, hi = matrix.getOwnershipRange()
    columns = np.arange(expected.shape[1], dtype=PETSc.IntType)
    rows = np.arange(lo, hi, dtype=PETSc.IntType)
    local = (matrix.getValues(rows, columns) if hi > lo
             else np.empty((0, len(columns)), dtype=complex))
    actual = np.concatenate(matrix.getComm().tompi4py().allgather(local))
    info = matrix.getInfo(PETSc.Mat.InfoType.GLOBAL_SUM)
    assert_allclose(actual, expected, rtol=0, atol=1e-14)
    assert info['mallocs'] == 0


@pytest.mark.parametrize('distributed', [False, True])
@pytest.mark.parametrize('size', [1, 16])
@pytest.mark.parametrize('hint', [None, 16])
def test_dense_one_electron_operator_preallocation(distributed, size, hint):
    """All-to-all hopping needs more than the default five remote entries/row.

    The size-one case also exercises empty ownership ranges under MPI.
    In the one-electron basis the many-body matrix equals the input hopping.
    """
    comm = PETSc.COMM_WORLD if distributed else PETSc.COMM_SELF
    rng = np.random.default_rng(367)
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    expected = raw + raw.conj().T
    matrix = build_op(
        expected, None, FockBasisSpec.from_args(size, 1), backend='petsc',
        use_numba=False, backend_kws={'comm': comm, 'nnz_guess_per_row': hint},
    )
    try:
        assert_matrix_values_and_allocations(matrix, expected)
    finally:
        matrix.destroy()


@pytest.mark.parametrize('distributed', [False, True])
@pytest.mark.parametrize('shape', [(24, 16), (1, 2), (2, 1)])
def test_rectangular_preallocation_with_column_owned_insertion(distributed, shape):
    """Transition-shaped matrices can have different row/column ownership."""
    comm = PETSc.COMM_WORLD if distributed else PETSc.COMM_SELF
    nl, nr = shape
    matrix = _create_matrix(PETSc, comm, nl, nr, nr)
    expected = np.arange(1, nl + 1)[:, None] + 1j * np.arange(1, nr + 1)[None, :]
    try:
        lo, hi = matrix.getOwnershipRangeColumn()
        rows = np.repeat(np.arange(nl), hi - lo)
        columns = np.tile(np.arange(lo, hi), nl)
        values = expected[rows, columns]
        _insert_entries(matrix, rows, columns, values, PETSc)
        _finalize_matrix(matrix, None)
        assert_matrix_values_and_allocations(matrix, expected)
    finally:
        matrix.destroy()
