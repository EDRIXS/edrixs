"""Backend implementation which delegates calculations to f2py solvers.

This code is currently a somewhat cumbersome bridge between the revised frontend
and the older Fortran. In the longer term, the Fortran will either be edited to
generate a cleaner interface, or the Fortran may be removed from the codebase.

The Fortran programs communicate exclusively through their conventional input
and output files in the current working directory.  ``get_ops`` writes the whole
problem there -- ``hopping_*.in``, ``coulomb_*.in``, ``fock_*.in``,
``config.in`` and ``transop_components.in`` -- and hands back plain marker
objects; :func:`ed`, :func:`xas` and :func:`rixs` operate on those files and
read ``eigvals.dat`` / ``*_poles.*`` back.
"""

from __future__ import annotations

from pathlib import Path
import traceback
import warnings

import numpy as np
import scipy.sparse as sp

from ..fock_basis import FockBasisSpec, build_fock_basis
from ..photon_transition import dipole_polvec_xas, dipole_polvec_rixs, quadrupole_polvec, unit_wavevector
from ..poles import get_spectra_from_poles, merge_pole_dicts
from .iostream_fortran import (
    read_poles_from_file, write_config, write_emat, write_umat,
)
from .options import OPTIONS, validate_options

__all__ = [
    'FortranDiskOperator', 'get_ops_disk', 'write_problem', 'ed_fortran',
    'xas_fortran', 'rixs_fortran', 'owns_operator_fortran',
]


class FortranDiskOperator:
    """Stateless marker identifying the disk-backed Fortran backend.

    The handle contains no operator data or path.  Subsequent calls using
    ``backend='fortran'`` read the native input and output files from the
    current working directory.
    """

    def __repr__(self):
        return f"<Fortran disk handle: using files in {Path.cwd()}>"

    __str__ = __repr__


def owns_operator_fortran(operator):
    return isinstance(operator, FortranDiskOperator)


def _read_config():
    """Recover ``(num_val_orbs, num_core_orbs)`` from the native ``config.in``.

    ``config.in`` is the single on-disk source of truth for the orbital counts,
    written by :func:`write_problem`; ``ed``/``xas``/``rixs`` read it back rather
    than carrying the numbers on the marker objects.
    """
    values = {}
    for line in Path('config.in').read_text().splitlines():
        key, sep, val = line.partition('=')
        if sep:
            values[key.strip()] = val.strip()
    return int(values['num_val_orbs']), int(values['num_core_orbs'])


def _communicator(options):
    """Return the communicator used by the disk-backed Fortran solvers."""
    if 'comm' in options:
        comm = options['comm']
    else:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    for method in ('Get_rank', 'Get_size', 'Barrier', 'bcast', 'py2f'):
        if not hasattr(comm, method):
            raise TypeError("backend_kws['comm'] must be an mpi4py communicator")
    return comm


def _root_collective(comm, action):
    """Run ``action`` on rank zero and return/broadcast its result everywhere."""
    if comm.Get_size() == 1:
        return action()

    comm.Barrier()
    payload = None
    if comm.Get_rank() == 0:
        try:
            payload = (True, action())
        except BaseException:
            payload = (False, traceback.format_exc())
    success, value = comm.bcast(payload, root=0)
    comm.Barrier()
    if not success:
        raise RuntimeError("Fortran backend rank-zero operation failed:\n{}".format(value))
    return value


def _basis_spec(basis, name):
    spec = basis if hasattr(basis, 'shapes') else getattr(basis, 'spec', None)
    if spec is None:
        raise ValueError(f"the Fortran backend needs a FockBasisSpec for {name}")
    return spec


def _write_fock_basis(basis, num_core_orbs, fname):
    """Mask off core orbitals and write unique, sorted native valence states.

    Python and Fortran both store orbital zero in the least significant bit.
    Core-hole configurations collapse onto the same intermediate valence state.
    """
    basis = build_fock_basis(basis)
    valence_mask = (1 << (basis.norbs - num_core_orbs)) - 1
    states = {basis.decode(i) & valence_mask for i in range(len(basis))}
    with Path(fname).open('w') as stream:
        stream.write(f'{len(states)}\n')
        for state in sorted(states):
            stream.write(f'{state}\n')


def _write_fock_sector(nv, nocc, nc, fname, *, core_hole=False):
    """Stage a legacy fixed-occupancy sector through the full-core convention."""
    basis = FockBasisSpec.from_args(nv, nocc, nc, nc - int(core_hole))
    _write_fock_basis(basis, nc, fname)


def _reduce_initial_integrals(emat, umat, num_val_orbs):
    """Contract a filled core for ``U[p,q,r,s] c†p c†q cr cs``.

    Return independent valence one-/two-body arrays and the core constant.
    Flattened sparse Coulomb tensors stay sparse throughout the reduction.
    """
    nv = num_val_orbs
    ntot = len(emat)
    h = np.array(emat[:nv, :nv], dtype=complex, copy=True)
    core_energy = np.trace(emat[nv:, nv:])
    if sp.issparse(umat):
        coo = umat.tocoo(copy=True)
        p, q = coo.row // ntot, coo.row % ntot
        r, s = coo.col // ntot, coo.col % ntot
        values = coo.data
        valence = (p < nv) & (q < nv) & (r < nv) & (s < nv)
        u = sp.coo_matrix(
            (values[valence], (p[valence] * nv + q[valence],
                               r[valence] * nv + s[valence])),
            shape=(nv * nv, nv * nv),
        ).tocsr()
        for mask, row, col, sign in (
            ((p < nv) & (q >= nv) & (q == r) & (s < nv), p, s, 1),
            ((p < nv) & (q >= nv) & (q == s) & (r < nv), p, r, -1),
            ((p >= nv) & (q < nv) & (p == r) & (s < nv), q, s, -1),
            ((p >= nv) & (q < nv) & (p == s) & (r < nv), q, r, 1),
        ):
            np.add.at(h, (row[mask], col[mask]), sign * values[mask])
        core = (p >= nv) & (q >= nv) & (r >= nv) & (s >= nv)
        core_energy += values[core & (p == s) & (q == r)].sum()
        core_energy -= values[core & (p == r) & (q == s)].sum()
    else:
        u = np.array(umat[:nv, :nv, :nv, :nv], copy=True)
        for a in range(nv, ntot):
            h += (umat[:nv, a, a, :nv] - umat[:nv, a, :nv, a]
                  - umat[a, :nv, a, :nv] + umat[a, :nv, :nv, a])
            for b in range(nv, ntot):
                core_energy += umat[a, b, b, a] - umat[a, b, a, b]
    return h, u, core_energy


def _write_initial_integrals(emat, umat, num_val_orbs, nocc, tol=1e-10,
                             *, include_core_energy=True):
    """Reduce and stage the initial integrals, without modifying caller data."""
    h, u, core_energy = _reduce_initial_integrals(emat, umat, num_val_orbs)
    if include_core_energy:
        if nocc == 0:
            if abs(core_energy) > tol:
                raise ValueError(
                    'native Fortran cannot represent nonzero core energy '
                    'in a zero-valence sector'
                )
        else:
            h += np.eye(num_val_orbs) * (core_energy / nocc)
    write_emat(h, 'hopping_i.in')
    write_umat(u, 'coulomb_i.in', tol)
    return core_energy


def _write_transition_components(trans_mat, fname):
    """Write dense transition components in the backend's private format."""
    trans_mat = np.asarray(trans_mat, dtype=complex)
    if trans_mat.ndim != 3:
        raise ValueError("trans_mat must be a three-dimensional array")

    with Path(fname).open('w') as stream:
        for component, row, col in np.ndindex(trans_mat.shape):
            value = trans_mat[component, row, col]
            stream.write(
                f"{component + 1:10d}    {row + 1:10d}    {col + 1:10d}    "
                f"{value.real:.15f}    {value.imag:.15f}    \n"
            )


def _make_disk_operators(num_trans_ops):
    """Return Fortran disk handles for a problem's staged operators."""
    return (
        FortranDiskOperator(),
        FortranDiskOperator(),
        [FortranDiskOperator() for _ in range(num_trans_ops)],
    )


def write_problem(emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat,
                  *, backend_kws=None):
    """Write the complete native Fortran problem input and return disk handles.

    Every native file is written to the current working directory.  Under an
    MPI communicator only rank zero writes; all ranks share one working
    directory, so no path needs to be exchanged.
    """
    options = validate_options('get_ops', backend_kws)
    comm = _communicator(options)
    spec_i = _basis_spec(basis_i, 'basis_i')
    num_core_orbs = spec_i.shapes[-1][0]
    num_val_orbs = spec_i.norbs - num_core_orbs
    v_noccu = sum(nocc for _, nocc in spec_i.shapes[:-1])
    trans_mat = np.asarray(trans_mat)
    if trans_mat.ndim != 3:
        raise ValueError("trans_mat must be a three-dimensional array")

    def write_inputs():
        write_emat(np.asarray(emat_n), 'hopping_n.in')
        tol = options.get('tol', OPTIONS['get_ops']['tol'].default)
        _write_initial_integrals(emat_i, umat_i, num_val_orbs, v_noccu, tol)
        write_umat(umat_n, 'coulomb_n.in', tol)
        _write_fock_basis(basis_i, num_core_orbs, 'fock_i.in')
        _write_fock_basis(basis_n, num_core_orbs, 'fock_n.in')
        _write_fock_basis(basis_i, num_core_orbs, 'fock_f.in')
        write_config(num_val_orbs=num_val_orbs, num_core_orbs=num_core_orbs)
        _write_transition_components(trans_mat, 'transop_components.in')

    _root_collective(comm, write_inputs)
    # The transition components stay on disk in transop_components.in; xas/rixs
    # reload them and write the polarization-specific native input files later.
    return _make_disk_operators(len(trans_mat))


def _solver(name):
    """Return an f2py solver lazily so importing EDRIXS stays lightweight."""
    from .. import fedrixs

    return getattr(fedrixs, name)


def _run_solver(solver, comm, outputs=()):
    """Call an f2py solver on every rank of the existing communicator."""
    outputs = tuple(Path(path) for path in outputs)

    def clear_outputs():
        for output in outputs:
            output.unlink(missing_ok=True)

    _root_collective(comm, clear_outputs)
    rank = comm.Get_rank()
    solver(comm.py2f(), rank, comm.Get_size())


def ed_fortran(hmat_i, num_evals=1, *, shift=0.0, backend_kws=None):
    if shift != 0:
        warnings.warn(
            'Fortran ED ignores shift; continuing without applying the shift.',
            UserWarning, stacklevel=2,
        )
    options = validate_options('ed', backend_kws)
    nvector = int(options.get('nvector', num_evals))
    if nvector > num_evals:
        raise ValueError("backend_kws['nvector'] cannot exceed num_evals")
    comm = _communicator(options)

    def prepare_ed():
        num_val_orbs, num_core_orbs = _read_config()
        write_config(
            ed_solver=options.get(
                'ed_solver', OPTIONS['ed']['ed_solver'].default
            ),
            num_val_orbs=num_val_orbs,
            num_core_orbs=num_core_orbs,
            neval=num_evals,
            nvector=nvector,
            ncv=options.get('ncv', max(3, num_evals + 2)),
            idump=options.get('idump', OPTIONS['ed']['idump'].default),
            maxiter=options.get('maxiter', OPTIONS['ed']['maxiter'].default),
            min_ndim=options.get(
                'min_ndim', OPTIONS['ed']['min_ndim'].default
            ),
            eigval_tol=options.get(
                'eigval_tol', OPTIONS['ed']['eigval_tol'].default
            ),
        )

    def read_eigenvalues():
        data = np.loadtxt('eigvals.dat', ndmin=2)
        if data.shape[0] < num_evals:
            raise RuntimeError("ed_fsolver produced fewer eigenvalues than requested")
        return np.asarray(data[:num_evals, 1], dtype=float)

    _root_collective(comm, prepare_ed)
    _run_solver(_solver('ed_fsolver'), comm, outputs=['eigvals.dat'])
    return _root_collective(comm, read_eigenvalues), FortranDiskOperator()


def _read_transition_components():
    """Read the polarization components from the current working directory."""
    num_val_orbs, num_core_orbs = _read_config()
    ntot = num_val_orbs + num_core_orbs
    rows = np.loadtxt('transop_components.in', ndmin=2)
    values = rows[:, 3] + 1j * rows[:, 4]
    npol = values.size // (ntot * ntot)
    return values.reshape(npol, ntot, ntot)


def get_ops_disk():
    """Return operator handles for the Fortran problem in this directory.

    The native ``config.in`` and ``transop_components.in`` files are read from
    the current working directory to recover the number of transition
    components.  No problem files are written or modified.
    """
    trans_mat = _read_transition_components()
    return _make_disk_operators(len(trans_mat))


def _gamma(value, mesh):
    return np.full(len(mesh), value, dtype=float) if np.isscalar(value) else np.asarray(value, dtype=float)


def xas_fortran(eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c=0.1,
                thin=1.0, phi=0.0, pol_type=None, temperature=1.0,
                scatter_axis=None, backend_kws=None):
    """Run XAS collectively on the existing MPI communicator."""
    options = validate_options('xas', backend_kws)
    num_gs = int(options.get('num_gs', len(eval_i)))
    if num_gs > len(eval_i):
        raise ValueError("backend_kws['num_gs'] cannot exceed len(eval_i)")
    nkryl = int(options.get('nkryl', OPTIONS['xas']['nkryl'].default))
    comm = _communicator(options)

    components = _root_collective(comm, _read_transition_components)
    if pol_type is None:
        pol_type = [('isotropic', 0)]
    scatter_axis = np.eye(3) if scatter_axis is None else np.asarray(scatter_axis)

    def prepare_xas():
        num_val_orbs, num_core_orbs = _read_config()
        write_config(
            num_val_orbs=num_val_orbs,
            num_core_orbs=num_core_orbs, num_gs=num_gs, nkryl=nkryl,
        )

    _root_collective(comm, prepare_xas)
    xas_fsolver = _solver('xas_fsolver')
    gamma = _gamma(gamma_c, ominc)
    result, poles = np.zeros((len(ominc), len(pol_type))), []
    for ipol, (kind, alpha) in enumerate(pol_type):
        if kind.strip() == 'isotropic':
            vectors = [np.eye(len(components))[i] for i in range(len(components))]
        else:
            vectors = [
                np.asarray(
                    dipole_polvec_xas(thin, phi, alpha, scatter_axis, kind),
                    complex,
                )
            ]
        item_poles = []
        for vector in vectors:
            if len(components) == 5 and kind.strip() != 'isotropic':
                vector = quadrupole_polvec(vector, unit_wavevector(thin, phi, scatter_axis, 'in'))
            transop = np.tensordot(vector, components, axes=(0, 0))
            _root_collective(
                comm,
                lambda transop=transop: write_emat(transop, 'transop_xas.in'),
            )
            pole_files = [f'xas_poles.{i + 1}' for i in range(num_gs)]
            _run_solver(xas_fsolver, comm, outputs=pole_files)
            pole = _root_collective(comm, lambda: read_poles_from_file(pole_files))
            item_poles.append(pole)
            result[:, ipol] += get_spectra_from_poles(pole, ominc, gamma, temperature)
        result[:, ipol] /= len(vectors)
        poles.append(merge_pole_dicts(item_poles) if len(item_poles) > 1 else item_poles[0])
    return result


def rixs_fortran(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
                 *, gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0,
                 phi=0.0, pol_type=None, temperature=1.0, scatter_axis=None,
                 skip_gs=False, return_poles=False, backend_kws=None):
    """Run RIXS collectively on the existing MPI communicator."""
    options = validate_options('rixs', backend_kws)
    num_gs = int(options.get('num_gs', len(eval_i)))
    if num_gs > len(eval_i):
        raise ValueError("backend_kws['num_gs'] cannot exceed len(eval_i)")
    nkryl = int(options.get('nkryl', OPTIONS['rixs']['nkryl'].default))
    comm = _communicator(options)

    components = _root_collective(comm, _read_transition_components)
    if skip_gs:
        raise NotImplementedError("skip_gs is not supported by the disk-backed Fortran solver")
    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    scatter_axis = np.eye(3) if scatter_axis is None else np.asarray(scatter_axis)
    gamma_in, gamma_out = _gamma(gamma_c, ominc), _gamma(gamma_f, eloss)
    rixs_fsolver = _solver('rixs_fsolver')
    result, poles = np.zeros((len(ominc), len(eloss), len(pol_type))), []
    for iom, omega in enumerate(ominc):
        def prepare_rixs():
            num_val_orbs, num_core_orbs = _read_config()
            write_config(
                num_val_orbs=num_val_orbs,
                num_core_orbs=num_core_orbs, num_gs=num_gs, nkryl=nkryl,
                linsys_max=options.get(
                    'linsys_max', OPTIONS['rixs']['linsys_max'].default
                ),
                linsys_tol=options.get(
                    'linsys_tol', OPTIONS['rixs']['linsys_tol'].default
                ), omega_in=omega,
                gamma_in=gamma_in[iom],
            )

        _root_collective(comm, prepare_rixs)
        row = []
        for ipol, (it, alpha, jt, beta) in enumerate(pol_type):
            vin, vout = dipole_polvec_rixs(
                thin, thout, phi, alpha, beta, scatter_axis, (it, jt),
            )
            if len(components) == 5:
                vin = quadrupole_polvec(vin, unit_wavevector(thin, phi, scatter_axis, 'in'))
                vout = quadrupole_polvec(vout, unit_wavevector(thout, phi, scatter_axis, 'out'))
            transop_in = np.tensordot(vin, components, axes=(0, 0))
            transop_out = np.conj(np.tensordot(vout, components, axes=(0, 0)).T)

            def write_transops():
                write_emat(transop_in, 'transop_rixs_i.in')
                write_emat(transop_out, 'transop_rixs_f.in')

            _root_collective(comm, write_transops)
            pole_files = [f'rixs_poles.{i + 1}' for i in range(num_gs)]
            _run_solver(rixs_fsolver, comm, outputs=pole_files)
            pole = _root_collective(comm, lambda: read_poles_from_file(pole_files))
            row.append(pole)
            result[iom, :, ipol] = get_spectra_from_poles(
                pole, eloss, gamma_out, temperature,
            )
        poles.append(row)
    return (result, poles) if return_poles else result
