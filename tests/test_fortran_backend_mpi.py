"""MPI integration coverage for the disk-backed Fortran backend."""

import os
import shutil
import subprocess
import sys

import pytest

from edrixs.fortran_backend import fortran_backend
from edrixs.models import model_1v1c


pytestmark = pytest.mark.integration


class _RootComm:
    """Minimal single-process stand-in for the rank-zero collective path."""

    def __init__(self):
        self.barriers = 0
        self.broadcasts = []

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 2

    def Barrier(self):
        self.barriers += 1

    def bcast(self, value, root):
        assert root == 0
        self.broadcasts.append(value)
        return value

    def py2f(self):
        return 42


class _NonrootComm:
    """Stand-in for a non-root rank receiving rank-zero broadcasts."""

    def __init__(self):
        self.barriers = 0
        self.broadcasts = 0

    def Get_rank(self):
        return 1

    def Get_size(self):
        return 2

    def Barrier(self):
        self.barriers += 1

    def bcast(self, value, root):
        assert root == 0
        self.broadcasts += 1
        # The only bcast reports the successful root-only write result.
        return (True, None)

    def py2f(self):
        return 43


def test_parent_mpi_root_collective_broadcasts_values_and_errors():
    comm = _RootComm()
    assert fortran_backend._root_collective(comm, lambda: 'root result') == 'root result'
    assert comm.barriers == 2
    assert comm.broadcasts == [(True, 'root result')]

    with pytest.raises(RuntimeError, match='rank-zero operation failed'):
        fortran_backend._root_collective(comm, lambda: (_ for _ in ()).throw(ValueError('boom')))


def test_parent_mpi_nonroot_does_not_write_problem_files(tmp_path, monkeypatch):
    """Only rank zero writes native inputs; non-root returns disk handles."""
    monkeypatch.chdir(tmp_path)
    comm = _NonrootComm()

    def forbidden_write(*args, **kwargs):
        raise AssertionError('non-root wrote a native input file')

    monkeypatch.setattr(fortran_backend, 'write_emat', forbidden_write)
    problem = model_1v1c(('s', 's'), v_noccu=1)
    hmat_i, hmat_n, transitions = fortran_backend.write_problem(
        *problem, backend_kws={'comm': comm},
    )
    assert fortran_backend.owns_operator_fortran(hmat_i)
    assert fortran_backend.owns_operator_fortran(hmat_n)
    assert len(transitions) == 5
    assert comm.barriers == 2


def test_run_solver_reuses_existing_communicator():
    comm = _RootComm()
    calls = []

    fortran_backend._run_solver(
        lambda fcomm, rank, size: calls.append((fcomm, rank, size)),
        comm,
    )

    assert calls == [(42, 0, 2)]


@pytest.mark.skipif(
    shutil.which('mpirun') is None,
    reason='requires mpirun',
)
def test_fortran_backend_reuses_existing_mpi_communicator(tmp_path):
    """All Python ranks call the f2py solvers on their shared communicator."""
    script = tmp_path / 'collective_fortran.py'
    script.write_text(
        """
import numpy as np
from mpi4py import MPI
from edrixs.models import model_1v1c
from edrixs.solvers import ed, get_ops, rixs, xas

comm = MPI.COMM_WORLD
problem = model_1v1c(('s', 's'), v_noccu=1)
hmat_i, hmat_n, transitions = get_ops(*problem, backend='fortran')
eval_i, evec_i = ed(hmat_i, num_evals=1)
absorption = xas(
    eval_i, evec_i, hmat_n, transitions, np.array([0.0]),
    backend_kws={'nkryl': 2},
)
scattering = rixs(
    eval_i, evec_i, hmat_i, hmat_n, transitions,
    np.array([0.0]), np.array([0.0]), backend_kws={'nkryl': 2},
)
for value in (eval_i, absorption, scattering):
    values = comm.allgather(value)
    assert all(np.array_equal(values[0], item, equal_nan=True) for item in values[1:])
"""
    )
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join(
        entry for entry in [env.get('PYTHONPATH', ''), *sys.path] if entry
    )
    subprocess.run(
        ['mpirun', '--oversubscribe', '-np', '2', sys.executable, str(script)],
        cwd=tmp_path,
        env=env,
        check=True,
    )
