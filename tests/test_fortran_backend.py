"""Tests for the disk-backed f2py Fortran backend."""

from pathlib import Path

import numpy as np
import scipy.sparse as sp

import edrixs
from edrixs.fortran_backend import fortran_backend
from edrixs.fortran_backend.isostream_fortran import write_config, write_umat
from edrixs.models import model_1v1c
from edrixs.solvers import ed, get_ops, rixs, xas


def _write_poles(stem):
    Path(f'{stem}.1').write_text(
        'npoles 1\n'
        'eigval 0.0\n'
        'norm 1.0\n'
        '1 1.0 0.0\n'
    )


def test_sparse_umat_is_written_without_densifying(tmp_path, monkeypatch):
    umat = np.zeros((2, 2, 2, 2), dtype=complex)
    umat[0, 1, 1, 0] = 1.25 - 0.5j
    umat[1, 0, 0, 1] = -0.75 + 0.25j
    sparse_umat = sp.csr_matrix(umat.reshape(4, 4))
    dense_file = tmp_path / 'dense.in'
    sparse_file = tmp_path / 'sparse.in'
    compatibility_file = tmp_path / 'compatibility.in'

    write_umat(umat, dense_file)
    edrixs.write_umat(umat, compatibility_file)

    def forbidden_toarray(*args, **kwargs):
        raise AssertionError('sparse Coulomb data was densified')

    monkeypatch.setattr(sp.csr_matrix, 'toarray', forbidden_toarray)
    write_umat(sparse_umat, sparse_file)

    assert sparse_file.read_text() == dense_file.read_text()
    assert compatibility_file.read_text() == dense_file.read_text()


def test_transition_components_round_trip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    components = np.arange(20).reshape(5, 2, 2).astype(complex)
    components += 0.5j * components
    write_config(num_val_orbs=1, num_core_orbs=1)

    fortran_backend._write_transition_components(
        components, 'transop_components.in',
    )

    np.testing.assert_allclose(
        fortran_backend._read_transition_components(), components,
    )


def test_fortran_backend_uses_native_files_and_f2py_solvers(tmp_path, monkeypatch):
    """The public API delegates all numerical work to the f2py solvers."""
    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_solver(name):
        def run(fcomm, rank, size):
            calls.append((name, fcomm, rank, size))
            if name == 'ed_fsolver':
                Path('eigvals.dat').write_text('1 -0.25\n')
            elif name == 'xas_fsolver':
                _write_poles('xas_poles')
            elif name == 'rixs_fsolver':
                _write_poles('rixs_poles')

        return run

    monkeypatch.setattr(
        'edrixs.fortran_backend.fortran_backend._solver', fake_solver,
    )
    problem = model_1v1c(('s', 's'), v_noccu=1)
    hmat_i, hmat_n, transitions = get_ops(*problem, backend='fortran')
    assert len(transitions) == 5
    handles = (hmat_i, hmat_n, *transitions)
    assert all(type(handle) is type(hmat_i) for handle in handles)
    assert all(not vars(handle) for handle in handles)

    assert str(tmp_path) in repr(hmat_i)
    for filename in (
        'hopping_i.in', 'hopping_n.in', 'coulomb_i.in', 'coulomb_n.in',
        'fock_i.in', 'fock_n.in', 'fock_f.in', 'config.in',
    ):
        assert (tmp_path / filename).is_file()

    eval_i, evec_i = ed(hmat_i)
    np.testing.assert_allclose(eval_i, [-0.25])
    assert type(evec_i) is type(hmat_i)
    assert not vars(evec_i)
    assert str(tmp_path) in repr(evec_i)

    absorption = xas(
        eval_i, evec_i, hmat_n, transitions, np.array([0.0]),
    )
    scattering = rixs(
        eval_i, evec_i, hmat_i, hmat_n, transitions,
        np.array([0.0]), np.array([0.0]),
    )

    assert absorption.shape == (1, 1)
    assert scattering.shape == (1, 1, 1)
    names = [name for name, fcomm, rank, size in calls]
    assert names == (
        ['ed_fsolver'] + ['xas_fsolver'] * len(transitions) + ['rixs_fsolver']
    )
    assert all(rank == 0 and size == 1 for name, fcomm, rank, size in calls)
