"""Readers and writers for the native EDRIXS Fortran file formats."""

from __future__ import annotations

from math import isqrt
from pathlib import Path

import numpy as np
import scipy.sparse as sp

__all__ = [
    'write_emat', 'write_umat', 'write_config', 'read_poles_from_file',
]


def write_emat(emat, fname, tol=1e-12, fmt_int='{:10d}',
               fmt_float='{:.15f}'):
    """Write a dense one-body matrix in the native sparse input format."""
    emat = np.asarray(emat)
    rows, cols = np.nonzero(np.abs(emat) > tol)
    entries = np.stack((rows, cols), axis=-1)
    space = '    '
    fmt_string = (fmt_int + space) * 2 + (fmt_float + space) * 2 + '\n'

    with Path(fname).open('w') as stream:
        if not len(entries):
            stream.write(f'{1:10d}\n')
            stream.write(fmt_string.format(1, 1, 0.0, 0.0))
            return

        stream.write(f'{len(entries):20d}\n')
        for row, col in entries:
            value = emat[row, col]
            stream.write(fmt_string.format(
                row + 1, col + 1, value.real, value.imag,
            ))


def write_umat(umat, fname, tol=1e-12, fmt_int='{:10d}',
               fmt_float='{:.15f}'):
    """Write dense or flattened sparse Coulomb data in native rank-4 format.

    Sparse input uses the ``(n*n, n*n)`` convention and is streamed directly
    from its stored entries without allocating a dense rank-4 tensor.
    """
    space = '    '
    fmt_string = (fmt_int + space) * 4 + (fmt_float + space) * 2 + '\n'

    if sp.issparse(umat):
        norbs = isqrt(umat.shape[0])
        expected = (norbs * norbs, norbs * norbs)
        if umat.shape != expected:
            raise ValueError(
                "sparse umat has shape {}, expected the flattened "
                "(n*n, n*n) convention".format(umat.shape)
            )

        umat = umat.tocsr(copy=True)
        umat.sum_duplicates()
        umat = umat.tocoo()
        keep = np.abs(umat.data) > tol
        rows = umat.row[keep]
        cols = umat.col[keep]
        values = umat.data[keep]
        entries = []
        for row, col, value in zip(rows, cols, values):
            lorb, korb = divmod(int(row), norbs)
            jorb, iorb = divmod(int(col), norbs)
            entries.append((lorb, korb, jorb, iorb, value))
    else:
        umat = np.asarray(umat)
        if umat.ndim != 4 or len(set(umat.shape)) != 1:
            raise ValueError("dense umat must have shape (n, n, n, n)")
        entries = [
            (lorb, korb, jorb, iorb, umat[lorb, korb, jorb, iorb])
            for lorb, korb, jorb, iorb
            in zip(*np.nonzero(np.abs(umat) > tol))
        ]

    with Path(fname).open('w') as stream:
        if not entries:
            stream.write(f'{1:10d}\n')
            stream.write(fmt_string.format(1, 1, 1, 1, 0.0, 0.0))
            return

        stream.write(f'{len(entries):20d}\n')
        for lorb, korb, jorb, iorb, value in entries:
            stream.write(fmt_string.format(
                lorb + 1, korb + 1, jorb + 1, iorb + 1,
                value.real, value.imag,
            ))


def write_config(
        directory='.', ed_solver=1, num_val_orbs=2, num_core_orbs=2,
        neval=1, nvector=1, ncv=1, idump=True, num_gs=1, maxiter=500,
        linsys_max=1000, min_ndim=1000, nkryl=500, eigval_tol=1e-8,
        linsys_tol=1e-10, omega_in=0.0, gamma_in=0.1):
    """Write the native Fortran solver namelist to ``config.in``."""
    dump_vector = '.true.' if idump else '.false.'
    config = [
        '&control',
        'ed_solver=' + str(ed_solver),
        'num_val_orbs=' + str(num_val_orbs),
        'num_core_orbs=' + str(num_core_orbs),
        'neval=' + str(neval),
        'nvector=' + str(nvector),
        'ncv=' + str(ncv),
        'idump=' + dump_vector,
        'num_gs=' + str(num_gs),
        'maxiter=' + str(maxiter),
        'linsys_max=' + str(linsys_max),
        'min_ndim=' + str(min_ndim),
        'nkryl=' + str(nkryl),
        'eigval_tol=' + str(eigval_tol),
        'linsys_tol=' + str(linsys_tol),
        'omega_in=' + str(omega_in),
        'gamma_in=' + str(gamma_in),
        '&end',
    ]
    Path(directory, 'config.in').write_text(
        ''.join(item + '\n' for item in config)
    )


def read_poles_from_file(file_list):
    """Read native XAS or RIXS pole files into a pole dictionary."""
    pole_dict = {
        'npoles': [],
        'eigval': [],
        'norm': [],
        'alpha': [],
        'beta': [],
    }
    for fname in file_list:
        with Path(fname).open() as stream:
            neff = int(stream.readline().strip().split()[1])
            pole_dict['npoles'].append(neff)
            pole_dict['eigval'].append(
                float(stream.readline().strip().split()[1])
            )
            pole_dict['norm'].append(
                float(stream.readline().strip().split()[1])
            )

            alpha, beta = [], []
            for _ in range(neff):
                line = stream.readline().strip().split()
                alpha.append(float(line[1]))
                beta.append(float(line[2]))
            pole_dict['alpha'].append(alpha)
            pole_dict['beta'].append(beta)

    return pole_dict
