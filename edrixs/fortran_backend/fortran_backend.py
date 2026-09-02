"""Backend implementation which delegates calculations to ``*.x`` programs.

The Fortran programs communicate exclusively through their conventional input
and output files in the current working directory.  ``get_ops`` writes the whole
problem there -- ``hopping_*.in``, ``coulomb_*.in``, ``fock_*.in``,
``config.in`` and ``transop_components.in`` -- and hands back plain marker
objects; :func:`ed`, :func:`xas` and :func:`rixs` operate on those files and
read ``eigvals.dat`` / ``*_poles.*`` back.
"""

from __future__ import annotations

from pathlib import Path
import shlex
import subprocess
import time
import traceback

import numpy as np
import scipy.sparse as sp

from ..fock_basis import write_fock_dec_by_N
from ..iostream import write_config, write_emat, write_umat, write_tensor, read_poles_from_file
from ..photon_transition import dipole_polvec_xas, dipole_polvec_rixs, quadrupole_polvec, unit_wavevector
from ..plot_spectrum import get_spectra_from_poles, merge_pole_dicts

__all__ = [
    'FortranDiskOperator', 'FortranEigenvectors', 'write_problem',
    'ed_fortran', 'xas_fortran', 'rixs_fortran', 'owns_operator_fortran',
]


class FortranDiskOperator:
    """Placeholder for a Fortran operator held in native files on disk.

    ``get_ops(..., backend='fortran')`` returns these in place of ``hmat_i``,
    ``hmat_n`` and ``trans_ops``.  They carry no array data -- everything lives
    in the working directory -- and only name the role each on-disk operator
    plays.
    """

    def __init__(self, role):
        self.role = role

    def __repr__(self):
        return f"<Fortran {self.role}: data on disk in {Path.cwd()}>"

    __str__ = __repr__


class FortranEigenvectors(FortranDiskOperator):
    """Placeholder for the ``eigvec.*`` files written by ``ed.x``."""

    def __init__(self):
        super().__init__('eigenvectors')


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


def _options(backend_kws):
    if backend_kws is None:
        return {}
    if not hasattr(backend_kws, 'get'):
        raise TypeError("backend_kws must be a mapping or None")
    return dict(backend_kws)


def _communicator(options):
    """Return the parent communicator used to coordinate disk-backed work."""
    if 'comm' in options:
        comm = options['comm']
    else:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    for method in ('Get_rank', 'Get_size', 'Barrier', 'bcast'):
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


def _write_umat(umat, fname, tol=1e-12):
    """Write the Coulomb tensor in native rank-4 format via :func:`write_umat`.

    ``models`` may hand back the interaction as a flattened ``(n*n, n*n)``
    sparse matrix; expand it back to a rank-4 tensor so the shared writer can
    emit the standard ``coulomb_*.in`` file.
    """
    if sp.issparse(umat):
        norbs = round(np.sqrt(umat.shape[0]))
        if umat.shape != (norbs * norbs, norbs * norbs):
            raise ValueError("sparse umat must use the flattened (n*n, n*n) convention")
        umat = umat.toarray().reshape(norbs, norbs, norbs, norbs)
    write_umat(np.asarray(umat), str(fname), tol=tol)


def write_problem(emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat,
                  *, backend_kws=None):
    """Write the complete native Fortran problem input and return disk handles.

    Every native file is written to the current working directory.  Under an
    MPI parent communicator only rank zero writes; all ranks share one working
    directory, so no path needs to be exchanged.
    """
    options = _options(backend_kws)
    comm = _communicator(options)
    # ``basis_i`` is the valence-only initial sector; ``basis_n`` appends the
    # core-hole sector.  The native fock files span the valence orbitals only,
    # with one extra valence electron in the intermediate state.
    num_val_orbs, v_noccu = _basis_spec(basis_i, 'basis_i').shapes[0]
    num_core_orbs = _basis_spec(basis_n, 'basis_n').shapes[-1][0]
    trans_mat = np.asarray(trans_mat)
    if trans_mat.ndim != 3:
        raise ValueError("trans_mat must be a three-dimensional array")

    def write_inputs():
        write_emat(np.asarray(emat_i), 'hopping_i.in')
        write_emat(np.asarray(emat_n), 'hopping_n.in')
        _write_umat(umat_i, 'coulomb_i.in', options.get('tol', 1e-12))
        _write_umat(umat_n, 'coulomb_n.in', options.get('tol', 1e-12))
        write_fock_dec_by_N(num_val_orbs, v_noccu, 'fock_i.in')
        write_fock_dec_by_N(num_val_orbs, v_noccu + 1, 'fock_n.in')
        write_fock_dec_by_N(num_val_orbs, v_noccu, 'fock_f.in')
        write_config(num_val_orbs=num_val_orbs, num_core_orbs=num_core_orbs)
        write_tensor(np.asarray(trans_mat, dtype=complex), 'transop_components.in')

    _root_collective(comm, write_inputs)
    # The transition components stay on disk in transop_components.in; xas/rixs
    # reload them and write the polarization-specific native input files later.
    handles = [FortranDiskOperator('transition') for _ in range(len(trans_mat))]
    return (
        FortranDiskOperator('initial Hamiltonian'),
        FortranDiskOperator('intermediate Hamiltonian'),
        handles,
    )


def _run(program, options, comm, outputs=()):
    outputs = tuple(Path(path) for path in outputs)
    for output in outputs:
        output.unlink(missing_ok=True)
    command = options.get('command', program)
    if comm.Get_size() > 1:
        if isinstance(command, str):
            command = shlex.split(command)
        else:
            command = list(command)
        if not command:
            raise ValueError("Fortran executable command must not be empty")
        if Path(command[0]).name in {'mpirun', 'mpiexec', 'srun'}:
            raise ValueError(
                "Do not use an MPI launcher in parent-MPI mode; the Fortran "
                "backend spawns the native solver with comm.Get_size() ranks"
            )
        _spawn_native(command, comm.Get_size())
    else:
        # A string is intentionally run through the shell so callers can supply an
        # MPI launcher, e.g. ``'mpirun -np 4 ed.x'``.  A sequence avoids a shell
        # and is preferred when no launcher syntax is needed.  Either way the
        # native program runs in the current working directory.
        subprocess.run(command, check=True, shell=isinstance(command, str))

    deadline = time.monotonic() + float(options.get('output_timeout', 60.0))
    while outputs and not all(output.is_file() for output in outputs):
        if time.monotonic() >= deadline:
            missing = ', '.join(str(path) for path in outputs if not path.is_file())
            raise RuntimeError("Fortran solver did not produce expected output: {}".format(missing))
        time.sleep(0.01)


def _spawn_native(command, num_procs):
    """Spawn one native MPI solver group and wait for its clean disconnect.

    Kept separate from :func:`_run` so the parent-MPI control flow can be
    tested without requiring an MPI runtime capable of dynamic spawning.
    """
    from mpi4py import MPI

    children = MPI.COMM_SELF.Spawn(command[0], args=command[1:], maxprocs=num_procs)
    children.Disconnect()


def _check_handle(handle):
    if not owns_operator_fortran(handle):
        raise TypeError("Fortran backend requires operators returned by get_ops(..., backend='fortran')")
    return handle


def ed_fortran(hmat_i, num_evals=1, *, backend_kws=None):
    _check_handle(hmat_i)
    options = _options(backend_kws)
    comm = _communicator(options)
    nvector = int(options.get('nvector', num_evals))

    def run_ed():
        num_val_orbs, num_core_orbs = _read_config()
        write_config(ed_solver=options.get('ed_solver', 1),
                     num_val_orbs=num_val_orbs, num_core_orbs=num_core_orbs,
                     neval=num_evals, nvector=nvector, ncv=options.get('ncv', max(3, num_evals + 2)),
                     idump=options.get('idump', True), maxiter=options.get('maxiter', 500),
                     min_ndim=options.get('min_ndim', 1000), eigval_tol=options.get('eigval_tol', 1e-8))
        _run(options.get('ed_executable', 'ed.x'), options, comm,
             outputs=['eigvals.dat'])
        data = np.loadtxt('eigvals.dat', ndmin=2)
        if data.shape[0] < num_evals:
            raise RuntimeError("ed.x produced fewer eigenvalues than requested")
        return np.asarray(data[:num_evals, 1], dtype=float)

    return _root_collective(comm, run_ed), FortranEigenvectors()


def _transitions(trans_op):
    handles = list(trans_op)
    if not handles:
        raise ValueError("at least one transition operator is required")
    for item in handles:
        _check_handle(item)
    num_val_orbs, num_core_orbs = _read_config()
    ntot = num_val_orbs + num_core_orbs
    rows = np.loadtxt('transop_components.in', ndmin=2)
    values = rows[:, 3] + 1j * rows[:, 4]
    npol = values.size // (ntot * ntot)
    return values.reshape(npol, ntot, ntot)


def _gamma(value, mesh):
    return np.full(len(mesh), value, dtype=float) if np.isscalar(value) else np.asarray(value, dtype=float)


def xas_fortran(eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c=0.1,
                thin=1.0, phi=0.0, pol_type=None, temperature=1.0,
                scatter_axis=None, backend_kws=None):
    """Run XAS collectively when called from an MPI parent communicator."""
    options = _options(backend_kws)
    comm = _communicator(options)
    return _root_collective(
        comm,
        lambda: _xas_fortran_root(
            eval_i, evec_i, hmat_n, trans_op, ominc, gamma_c=gamma_c,
            thin=thin, phi=phi, pol_type=pol_type, temperature=temperature,
            scatter_axis=scatter_axis, options=options, comm=comm,
        ),
    )


def _xas_fortran_root(eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c,
                      thin, phi, pol_type, temperature, scatter_axis, options, comm):
    components = _transitions(trans_op)
    _check_handle(hmat_n)
    if pol_type is None:
        pol_type = [('isotropic', 0)]
    scatter_axis = np.eye(3) if scatter_axis is None else np.asarray(scatter_axis)
    num_gs = int(options.get('num_gs', len(eval_i)))
    nkryl = int(options.get('nkryl', 200))
    num_val_orbs, num_core_orbs = _read_config()
    write_config(
        num_val_orbs=num_val_orbs,
        num_core_orbs=num_core_orbs, num_gs=num_gs, nkryl=nkryl,
    )
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
            write_emat(np.tensordot(vector, components, axes=(0, 0)), 'transop_xas.in')
            pole_files = [f'xas_poles.{i + 1}' for i in range(num_gs)]
            _run(options.get('xas_executable', 'xas.x'), options, comm,
                 outputs=pole_files)
            pole = read_poles_from_file(pole_files)
            item_poles.append(pole)
            result[:, ipol] += get_spectra_from_poles(pole, ominc, gamma, temperature)
        result[:, ipol] /= len(vectors)
        poles.append(merge_pole_dicts(item_poles) if len(item_poles) > 1 else item_poles[0])
    return result


def rixs_fortran(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
                 *, gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0,
                 phi=0.0, pol_type=None, temperature=1.0, scatter_axis=None,
                 skip_gs=False, return_poles=False, backend_kws=None):
    """Run RIXS collectively when called from an MPI parent communicator."""
    options = _options(backend_kws)
    comm = _communicator(options)
    return _root_collective(
        comm,
        lambda: _rixs_fortran_root(
            eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
            gamma_c=gamma_c, gamma_f=gamma_f, thin=thin, thout=thout,
            phi=phi, pol_type=pol_type, temperature=temperature,
            scatter_axis=scatter_axis, skip_gs=skip_gs,
            return_poles=return_poles, options=options, comm=comm,
        ),
    )


def _rixs_fortran_root(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
                       *, gamma_c, gamma_f, thin, thout, phi, pol_type,
                       temperature, scatter_axis, skip_gs, return_poles,
                       options, comm):
    components = _transitions(trans_op)
    _check_handle(hmat_i)
    _check_handle(hmat_n)
    if skip_gs:
        raise NotImplementedError("skip_gs is not supported by the disk-backed Fortran solver")
    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    scatter_axis = np.eye(3) if scatter_axis is None else np.asarray(scatter_axis)
    num_gs = int(options.get('num_gs', len(eval_i)))
    nkryl = int(options.get('nkryl', 200))
    num_val_orbs, num_core_orbs = _read_config()
    gamma_in, gamma_out = _gamma(gamma_c, ominc), _gamma(gamma_f, eloss)
    result, poles = np.zeros((len(ominc), len(eloss), len(pol_type))), []
    for iom, omega in enumerate(ominc):
        write_config(
            num_val_orbs=num_val_orbs,
            num_core_orbs=num_core_orbs, num_gs=num_gs, nkryl=nkryl,
            linsys_max=options.get('linsys_max', 1000),
            linsys_tol=options.get('linsys_tol', 1e-10), omega_in=omega,
            gamma_in=gamma_in[iom],
        )
        row = []
        for ipol, (it, alpha, jt, beta) in enumerate(pol_type):
            vin, vout = dipole_polvec_rixs(thin, thout, phi, alpha, beta, scatter_axis, (it, jt))
            if len(components) == 5:
                vin = quadrupole_polvec(vin, unit_wavevector(thin, phi, scatter_axis, 'in'))
                vout = quadrupole_polvec(vout, unit_wavevector(thout, phi, scatter_axis, 'out'))
            write_emat(
                np.tensordot(vin, components, axes=(0, 0)),
                'transop_rixs_i.in',
            )
            write_emat(
                np.conj(np.tensordot(vout, components, axes=(0, 0)).T),
                'transop_rixs_f.in',
            )
            pole_files = [f'rixs_poles.{i + 1}' for i in range(num_gs)]
            _run(options.get('rixs_executable', 'rixs.x'), options, comm,
                 outputs=pole_files)
            pole = read_poles_from_file(pole_files)
            row.append(pole)
            result[iom, :, ipol] = get_spectra_from_poles(pole, eloss, gamma_out, temperature)
        poles.append(row)
    return (result, poles) if return_poles else result
