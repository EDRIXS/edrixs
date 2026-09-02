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


class _NonrootComm:
    """Stand-in for a non-root rank receiving rank-zero broadcasts."""

    def __init__(self, directory):
        self.directory = str(directory)
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
        # The first bcast supplies rank zero's working directory, and the
        # second supplies the successful root-only write result.
        return self.directory if self.broadcasts == 1 else (True, None)


def test_parent_mpi_root_collective_broadcasts_values_and_errors():
    comm = _RootComm()
    assert fortran_backend._root_collective(comm, lambda: 'root result') == 'root result'
    assert comm.barriers == 2
    assert comm.broadcasts == [(True, 'root result')]

    with pytest.raises(RuntimeError, match='rank-zero operation failed'):
        fortran_backend._root_collective(comm, lambda: (_ for _ in ()).throw(ValueError('boom')))


def test_parent_mpi_nonroot_does_not_write_problem_files(tmp_path, monkeypatch):
    """Only rank zero writes native inputs; non-root returns disk handles."""
    comm = _NonrootComm(tmp_path)

    def forbidden_write(*args, **kwargs):
        raise AssertionError('non-root wrote a native input file')

    monkeypatch.setattr(fortran_backend, 'write_emat', forbidden_write)
    problem = model_1v1c(('s', 's'), v_noccu=1)
    hmat_i, hmat_n, transitions = fortran_backend.write_problem(
        *problem, backend_kws={'comm': comm},
    )
    assert hmat_i.directory == tmp_path
    assert hmat_n.directory == tmp_path
    assert len(transitions) == 5
    assert comm.barriers == 2


def test_parent_mpi_run_uses_one_spawn_and_rejects_launchers(tmp_path, monkeypatch):
    comm = _RootComm()
    launches = []

    monkeypatch.setattr(
        fortran_backend,
        '_spawn_native',
        lambda command, num_procs: launches.append((command, num_procs)),
    )
    fortran_backend._run('fake-solver --flag', tmp_path, {}, comm)
    assert launches == [(['fake-solver', '--flag'], 2)]
    with pytest.raises(ValueError, match='Do not use an MPI launcher'):
        fortran_backend._run('mpirun -np 2 fake-solver', tmp_path, {}, comm)


@pytest.mark.skipif(
    not all(shutil.which(command) for command in ('mpirun', 'ed.x', 'xas.x', 'rixs.x')),
    reason='requires mpirun and the standalone Fortran solver commands',
)
def test_fortran_backend_collectively_spawns_native_solvers(tmp_path):
    """Two parent ranks launch one two-rank child group per native call."""
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
