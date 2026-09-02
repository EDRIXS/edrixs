"""Tests for the disk-backed standalone-Fortran backend."""

from pathlib import Path

import numpy as np

from edrixs.models import model_1v1c
from edrixs.solvers import ed, get_ops, rixs, xas


def _write_poles(stem):
    Path(f'{stem}.1').write_text(
        'npoles 1\n'
        'eigval 0.0\n'
        'norm 1.0\n'
        '1 1.0 0.0\n'
    )


def test_fortran_backend_uses_native_files_and_executables(tmp_path, monkeypatch):
    """The public API delegates all numerical work to the three commands."""
    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_run(command, check, shell):
        calls.append(command)
        program = command
        if program == 'fake-ed':
            Path('eigvals.dat').write_text('1 -0.25\n')
        elif program == 'fake-xas':
            _write_poles('xas_poles')
        elif program == 'fake-rixs':
            _write_poles('rixs_poles')

    monkeypatch.setattr('edrixs.fortran_backend.fortran_backend.subprocess.run', fake_run)
    problem = model_1v1c(('s', 's'), v_noccu=1)
    hmat_i, hmat_n, transitions = get_ops(*problem, backend='fortran')
    assert len(transitions) == 5

    assert 'on disk' in repr(hmat_i)
    for filename in (
        'hopping_i.in', 'hopping_n.in', 'coulomb_i.in', 'coulomb_n.in',
        'fock_i.in', 'fock_n.in', 'fock_f.in', 'config.in',
    ):
        assert (tmp_path / filename).is_file()

    eval_i, evec_i = ed(
        hmat_i, backend_kws={'ed_executable': 'fake-ed'}
    )
    np.testing.assert_allclose(eval_i, [-0.25])
    assert 'on disk' in repr(evec_i)

    absorption = xas(
        eval_i, evec_i, hmat_n, transitions, np.array([0.0]),
        backend_kws={'xas_executable': 'fake-xas'},
    )
    scattering = rixs(
        eval_i, evec_i, hmat_i, hmat_n, transitions,
        np.array([0.0]), np.array([0.0]),
        backend_kws={'rixs_executable': 'fake-rixs'},
    )

    assert absorption.shape == (1, 1)
    assert scattering.shape == (1, 1, 1)
    assert calls == (
        ['fake-ed'] + ['fake-xas'] * len(transitions) + ['fake-rixs']
    )
