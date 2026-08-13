"""PETSc backend interface.

This module implements the EDRIXS backend contract on top of petsc4py and
slepc4py. Importing EDRIXS does not require petsc4py; the heavy modules are
imported lazily inside the functions that need them, so simply importing this
module stays cheap and side-effect free.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

__all__ = [
    'owns_operator_petsc',
    'build_op_petsc',
    'ed_petsc', 'xas_petsc', 'rixs_petsc',
]


def _petsc_module():
    """
    Import petsc4py lazily and return its PETSc module.
    """
    try:
        from petsc4py import PETSc
    except ImportError as exc:
        raise ImportError(
            "The PETSc backend requires petsc4py and a working PETSc installation"
        ) from exc
    return PETSc


def _slepc_module():
    """
    Import slepc4py lazily and return its SLEPc module.
    """
    try:
        from slepc4py import SLEPc
    except ImportError as exc:
        raise ImportError(
            "The PETSc backend eigensolver requires slepc4py and a working "
            "SLEPc installation"
        ) from exc
    return SLEPc


def _backend_kws(backend_kws):
    """
    Validate and copy PETSc backend keyword arguments.
    """
    if backend_kws is None:
        return {}
    if not isinstance(backend_kws, Mapping):
        raise TypeError("backend_kws must be a mapping or None")
    return dict(backend_kws)


def owns_operator_petsc(operator):
    """
    Return whether ``operator`` is a petsc4py matrix or linear operator.
    """
    try:
        PETSc = _petsc_module()
    except ImportError:
        return False
    return isinstance(operator, (PETSc.Mat, PETSc.Vec))


def _not_implemented(operation):
    """
    Raise the standard PETSc-backend stub error.
    """
    raise NotImplementedError(
        "The PETSc backend contract is present, but {} has not yet been "
        "implemented".format(operation)
    )


# -----------------------------------------------------------------------------
# Many-body operator construction
# -----------------------------------------------------------------------------


def build_op_petsc(emat, umat, lb, rb=None, *, backend_kws=None):
    """Build a distributed PETSc many-body operator.

    Assemble ``H = sum_ij emat_ij f_i^dagger f_j
    + sum_lkji umat_lkji f_l^dagger f_k^dagger f_j f_i`` in the Fock basis
    ``lb`` (rows) and ``rb`` (columns). For a square operator omit ``rb``;
    ``lb`` is then used on both sides. For a rectangular transition operator
    ``rb`` is the column (right) basis and ``lb`` is the row (left) basis.

    Parameters
    ----------
    emat : 2d complex array or None
        One-body coefficients. ``None`` skips the one-body term.
    umat : 4d complex array or None
        Two-body Coulomb tensor. ``None`` skips the two-body term.
    lb : FockBinByN
        Left (row) many-body basis.
    rb : FockBinByN, optional
        Right (column) many-body basis. Defaults to ``lb``.
    backend_kws : mapping, optional
        Extra options. Recognized keys:

        - ``comm`` : PETSc communicator (default ``PETSc.COMM_WORLD``).
        - ``tol_e`` : threshold for retaining ``emat`` entries.
        - ``tol_u`` : threshold for retaining ``umat`` entries.
        - ``nnz_guess_per_row`` : preallocation hint.

    Returns
    -------
    petsc4py.PETSc.Mat
        The assembled many-body operator.
    """
    PETSc = _petsc_module()
    from .hash_basis_methods import get_H

    kws = _backend_kws(backend_kws)
    comm = kws.pop('comm', PETSc.COMM_WORLD)
    get_H_kws = {}
    for key in ('tol_e', 'tol_u', 'nnz_guess_per_row'):
        if key in kws:
            get_H_kws[key] = kws.pop(key)
    if kws:
        raise TypeError(
            "Unknown PETSc operator-construction options: {}".format(sorted(kws))
        )

    return get_H(comm, emat, umat, lb, rb=rb, **get_H_kws)


# -----------------------------------------------------------------------------
# Exact diagonalization (SLEPc)
# -----------------------------------------------------------------------------


def ed_petsc(hmat_i, num_evals=1, *, backend_kws=None):
    """Obtain the lowest-lying eigenpairs with SLEPc/PETSc.

    Diagonalize the Hermitian operator ``hmat_i`` for its ``num_evals``
    smallest-real eigenpairs using a Krylov-Schur solver.

    Parameters
    ----------
    hmat_i : petsc4py.PETSc.Mat
        Square Hermitian many-body Hamiltonian, e.g. from
        :func:`build_op_petsc`.
    num_evals : int, optional
        Number of lowest eigenpairs to return.
    backend_kws : mapping, optional
        Extra options. Recognized keys:

        - ``eigval_tol`` : convergence tolerance (default ``1e-8``).
        - ``maxiter`` : maximum solver iterations (default ``1000``).
        - ``ncv`` : number of column vectors (Krylov subspace size).
        - ``verbose`` : print convergence diagnostics (default ``False``).

    Returns
    -------
    eval_i : 1d float ndarray
        The ``num_evals`` lowest eigenvalues.
    evec_i : list of petsc4py.PETSc.Vec
        The corresponding eigenvectors (independent copies).
    """
    SLEPc = _slepc_module()

    kws = _backend_kws(backend_kws)
    eigval_tol = kws.pop('eigval_tol', 1e-8)
    maxiter = kws.pop('maxiter', 1000)
    ncv = kws.pop('ncv', None)
    verbose = kws.pop('verbose', False)
    if kws:
        raise TypeError("Unknown PETSc ED options: {}".format(sorted(kws)))

    num_evals = int(num_evals)
    if num_evals < 1:
        raise ValueError("num_evals must be a positive integer")

    comm = hmat_i.getComm()
    rank = comm.getRank()

    eps = SLEPc.EPS().create(comm=comm)
    eps.setOperators(hmat_i)
    eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
    eps.setProblemType(SLEPc.EPS.ProblemType.HEP)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
    eps.setTolerances(tol=eigval_tol, max_it=maxiter)
    if ncv is not None:
        eps.setDimensions(num_evals, int(ncv))
    else:
        eps.setDimensions(num_evals)

    eps.solve()

    nconv = eps.getConverged()
    if verbose and rank == 0:
        print("edrixs >>> PETSc ED converged {} of {} requested "
              "(reason={})".format(nconv, num_evals, eps.getConvergedReason()),
              flush=True)
    if nconv < num_evals:
        raise RuntimeError(
            "SLEPc converged only {} of the requested {} eigenpairs".format(
                nconv, num_evals
            )
        )

    eval_i = np.zeros(num_evals, dtype=float)
    evec_i = [None] * num_evals
    errors = np.zeros(num_evals, dtype=float)
    vr = hmat_i.getVecLeft()
    for index in range(num_evals):
        eigenvalue = eps.getEigenpair(index, vr)
        eval_i[index] = eigenvalue.real
        evec_i[index] = vr.copy()  # copy is crucial; vr is reused each pass
        errors[index] = eps.computeError(index, SLEPc.EPS.ErrorType.RELATIVE)

    if np.any(errors > eigval_tol):
        raise RuntimeError(
            "PETSc ED residuals exceed tolerance {}: {}".format(eigval_tol, errors)
        )

    return eval_i, evec_i


# -----------------------------------------------------------------------------
# Spectroscopy helpers
# -----------------------------------------------------------------------------


def _rixs_polarization_vectors(
        ntrans, thin, thout, phi, incoming_kind, alpha,
        outgoing_kind, beta, scatter_axis):
    """
    Return incoming and outgoing polarization vectors in transition-operator space.
    """
    from ..photon_transition import (
        dipole_polvec_rixs, quadrupole_polvec, unit_wavevector,
    )

    incoming, outgoing = dipole_polvec_rixs(
        thin, thout, phi, alpha, beta, scatter_axis,
        (incoming_kind, outgoing_kind),
    )
    if incoming_kind.lower() == 'isotropic':
        incoming = np.ones(3, dtype=complex) / np.sqrt(3.0)
    if outgoing_kind.lower() == 'isotropic':
        outgoing = np.ones(3, dtype=complex) / np.sqrt(3.0)

    incoming_vector = np.zeros(ntrans, dtype=complex)
    outgoing_vector = np.zeros(ntrans, dtype=complex)
    if ntrans == 3:
        incoming_vector[:] = incoming
        outgoing_vector[:] = outgoing
    elif ntrans == 5:
        incoming_wavevector = unit_wavevector(
            thin, phi, scatter_axis, direction='in'
        )
        outgoing_wavevector = unit_wavevector(
            thout, phi, scatter_axis, direction='out'
        )
        incoming_vector[:] = quadrupole_polvec(incoming, incoming_wavevector)
        outgoing_vector[:] = quadrupole_polvec(outgoing, outgoing_wavevector)
    else:
        raise ValueError("ntrans must be 3 or 5")
    return incoming_vector, outgoing_vector


def _apply_trans_combination(operators, coefficients, vector):
    """
    Apply ``sum_k coeff[k] * operators[k] @ vector`` (initial -> intermediate).
    """
    result = None
    for coefficient, operator in zip(coefficients, operators):
        if coefficient == 0:
            continue
        term = operator.createVecLeft()
        operator.mult(vector, term)
        term.scale(coefficient)
        if result is None:
            result = term
        else:
            result.axpy(1.0, term)
    if result is None:
        result = operators[0].createVecLeft()
        result.zeroEntries()
    return result


def _apply_trans_combination_adjoint(operators, coefficients, vector):
    """
    Apply ``sum_k coeff[k] * operators[k]^H @ vector`` (intermediate -> initial).
    """
    result = None
    for coefficient, operator in zip(coefficients, operators):
        if coefficient == 0:
            continue
        term = operator.createVecRight()
        operator.multHermitian(vector, term)
        term.scale(coefficient)
        if result is None:
            result = term
        else:
            result.axpy(1.0, term)
    if result is None:
        result = operators[0].createVecRight()
        result.zeroEntries()
    return result


def _xas_pole_dict(hmat_n, seeds, eval_i, nkryl, lanczos_tridiagonal):
    """
    Build an EDRIXS-compatible XAS pole dictionary from seed vectors.

    ``seeds[k]`` is the transition-operator-applied initial state associated
    with energy ``eval_i[k]`` and living in the intermediate Hilbert space.
    """
    pole_dict = {'npoles': [], 'eigval': [], 'norm': [],
                 'alpha': [], 'beta': []}
    for energy, seed in zip(eval_i, seeds):
        if seed.norm() == 0:
            alpha_i = np.array([0.0], dtype=float)
            beta_i = np.array([], dtype=float)
            norm_i = 0.0
        else:
            alpha_i, beta_i, norm_i = lanczos_tridiagonal(
                hmat_n, seed, nkryl=nkryl
            )
        pole_dict['npoles'].append(len(alpha_i))
        pole_dict['eigval'].append(float(energy))
        pole_dict['norm'].append(norm_i)
        pole_dict['alpha'].append(alpha_i)
        pole_dict['beta'].append(beta_i)
    return pole_dict


def xas_petsc(eval_i, evec_i, hmat_n, trans_op, ominc, *,
              gamma_c=0.1, thin=1.0, phi=0.0, pol_type=None,
              temperature=1.0, scatter_axis=None, backend_kws=None):
    """Calculate XAS with a PETSc Lanczos continued-fraction solver.

    ``eval_i`` and ``evec_i`` are the retained initial states, normally
    returned by :func:`ed_petsc`. For every retained initial state the
    transition operator is applied in the original Fock basis, and a Lanczos
    tridiagonalization of the intermediate Hamiltonian ``hmat_n`` generates
    the pole representation consumed by
    :func:`edrixs.plot_spectrum.get_spectra_from_poles`.

    Parameters
    ----------
    eval_i : 1d array
        Energies of the retained initial states.
    evec_i : sequence of petsc4py.PETSc.Vec
        Retained initial-state eigenvectors (e.g. from :func:`ed_petsc`).
        ``evec_i[k]`` corresponds to ``eval_i[k]``.
    hmat_n : petsc4py.PETSc.Mat
        Intermediate-state Hamiltonian.
    trans_op : sequence of petsc4py.PETSc.Mat
        Transition operators mapping the initial Hilbert space to the
        intermediate Hilbert space. The sequence length must be 3 (dipole)
        or 5 (quadrupole).
    ominc : 1d array
        Incident photon-energy grid.
    gamma_c : float or 1d array, optional
        Core-hole lifetime broadening. Scalar or same shape as ``ominc``.
    thin, phi : float, optional
        Incoming and azimuthal angles, in radians.
    pol_type : sequence of tuple, optional
        Incoming polarizations. Each entry is ``(kind, alpha)`` where ``kind``
        is ``'linear'``, ``'left'``, ``'right'``, or ``'isotropic'``. The
        default is a single isotropic polarization.
    temperature : float, optional
        Temperature in kelvin used for the Boltzmann weights of ``eval_i``.
    scatter_axis : (3, 3) array, optional
        Scattering coordinate axes. Defaults to the identity matrix.
    backend_kws : mapping, optional
        Extra options. Recognized keys:

        - ``nkryl`` : maximum intermediate-state Lanczos dimension
          (default ``200``).

    Returns
    -------
    xas : 2d ndarray
        Spectrum with shape ``(len(ominc), len(pol_type))``.
    poles : list of dict
        Per-polarization pole dictionaries.
    """
    _petsc_module()
    from .lanczos import lanczos_tridiagonal
    from ..plot_spectrum import get_spectra_from_poles, merge_pole_dicts
    from .._solvers_helpers import _expand_broadening
    from ..photon_transition import (
        dipole_polvec_xas, quadrupole_polvec, unit_wavevector,
    )

    kws = _backend_kws(backend_kws)
    nkryl = int(kws.pop('nkryl', 200))
    if kws:
        raise TypeError("Unknown PETSc XAS options: {}".format(sorted(kws)))

    eval_i = np.asarray(eval_i, dtype=float)
    ominc = np.asarray(ominc, dtype=float)
    evec_i = list(evec_i)
    num_gs = len(evec_i)
    if len(eval_i) != num_gs:
        raise ValueError("len(eval_i) must equal len(evec_i)")

    ntrans = len(trans_op)
    if ntrans not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )

    if pol_type is None:
        pol_type = [('isotropic', 0)]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.asarray(scatter_axis, dtype=float)
    if scatter_axis.shape != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

    n_om = len(ominc)
    gamma_core = _expand_broadening(gamma_c, n_om, 'gamma_c')
    kvec = unit_wavevector(thin, phi, scatter_axis, direction='in')

    xas = np.zeros((n_om, len(pol_type)), dtype=float)
    poles = []
    for it, (pt, alpha) in enumerate(pol_type):
        kind = pt.strip().lower()
        if kind == 'isotropic':
            # Sum over the Cartesian (or spherical) components individually.
            pole_dicts = []
            for k in range(ntrans):
                coeffs = np.zeros(ntrans, dtype=complex)
                coeffs[k] = 1.0
                seeds = [
                    _apply_trans_combination(trans_op, coeffs, evec_i[ig])
                    for ig in range(num_gs)
                ]
                pole_dict = _xas_pole_dict(
                    hmat_n, seeds, eval_i, nkryl, lanczos_tridiagonal
                )
                xas[:, it] += get_spectra_from_poles(
                    pole_dict, ominc, gamma_core, temperature
                ) / ntrans
                pole_dicts.append(pole_dict)
            poles.append(merge_pole_dicts(pole_dicts))
        elif kind in ('linear', 'left', 'right'):
            pol = dipole_polvec_xas(thin, phi, alpha, scatter_axis, pt)
            polvec = np.zeros(ntrans, dtype=complex)
            if ntrans == 3:
                polvec[:] = pol
            else:
                polvec[:] = quadrupole_polvec(pol, kvec)
            seeds = [
                _apply_trans_combination(trans_op, polvec, evec_i[ig])
                for ig in range(num_gs)
            ]
            pole_dict = _xas_pole_dict(
                hmat_n, seeds, eval_i, nkryl, lanczos_tridiagonal
            )
            poles.append(pole_dict)
            xas[:, it] = get_spectra_from_poles(
                pole_dict, ominc, gamma_core, temperature
            )
        else:
            raise ValueError("Unknown XAS polarization type: {}".format(pt))

    return xas, poles


def rixs_petsc(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
               gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
               pol_type=None, temperature=1.0, scatter_axis=None,
               skip_gs=False, return_poles=False, backend_kws=None):
    """Calculate RIXS with the PETSc Krylov correction-vector solver.

    For every incident energy, polarization, and retained initial state, the
    intermediate correction vector is obtained by solving the shifted linear
    system ``(omega + E_i + i*gamma_c - H_n) x = D_k |i>`` with a PETSc KSP
    (GMRES), and a Lanczos tridiagonalization of ``hmat_i`` seeded by
    ``D_k'^dagger x`` produces the pole representation consumed by
    :func:`edrixs.plot_spectrum.get_spectra_from_poles`.

    Parameters
    ----------
    eval_i : 1d array
        Energies of the retained initial states.
    evec_i : sequence of petsc4py.PETSc.Vec
        Retained initial-state eigenvectors (e.g. from :func:`ed_petsc`).
        ``evec_i[k]`` corresponds to ``eval_i[k]``.
    hmat_i, hmat_n : petsc4py.PETSc.Mat
        Initial/final and intermediate Hamiltonians.
    trans_op : sequence of petsc4py.PETSc.Mat
        Transition operators mapping the initial Hilbert space to the
        intermediate Hilbert space. The sequence length must be 3 (dipole)
        or 5 (quadrupole).
    ominc, eloss : 1d arrays
        Incident-energy and energy-loss grids.
    gamma_c, gamma_f : float or 1d array, optional
        Core-hole and final-state broadenings.
    thin, thout, phi : float, optional
        Scattering angles in radians.
    pol_type : sequence of tuple, optional
        Incoming and outgoing polarization specifications.
    temperature : float, optional
        Temperature in kelvin used for Boltzmann weights.
    scatter_axis : (3, 3) array, optional
        Scattering coordinate axes.
    skip_gs : bool, optional
        If true, project the retained initial-state subspace out of the
        final-state seed vector (removes the elastic line).
    return_poles : bool, optional
        Return the nested pole dictionaries together with the spectrum.
    backend_kws : mapping, optional
        Extra options. Recognized keys:

        - ``nkryl`` : maximum final-state Lanczos dimension (default ``200``).
        - ``linsys_tol`` : KSP absolute tolerance (default ``1e-10``).
        - ``linsys_max`` : maximum KSP iterations (default ``1000``).
        - ``ksp_type`` : KSP method (default ``'gmres'``).

    Returns
    -------
    rixs : 3d ndarray
        Spectrum with shape ``(len(ominc), len(eloss), len(pol_type))``.
    poles : list, optional
        Returned only when ``return_poles`` is true.
    """
    PETSc = _petsc_module()
    from .lanczos import lanczos_tridiagonal
    from ..plot_spectrum import get_spectra_from_poles
    from .._solvers_helpers import _expand_broadening

    kws = _backend_kws(backend_kws)
    nkryl = int(kws.pop('nkryl', 200))
    linsys_tol = float(kws.pop('linsys_tol', 1e-10))
    linsys_max = int(kws.pop('linsys_max', 1000))
    ksp_type = kws.pop('ksp_type', 'gmres')
    if kws:
        raise TypeError("Unknown PETSc RIXS options: {}".format(sorted(kws)))

    eval_i = np.asarray(eval_i, dtype=float)
    ominc = np.asarray(ominc, dtype=float)
    eloss = np.asarray(eloss, dtype=float)
    evec_i = list(evec_i)
    num_gs = len(evec_i)
    if len(eval_i) != num_gs:
        raise ValueError("len(eval_i) must equal len(evec_i)")

    ntrans = len(trans_op)
    if ntrans not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )

    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.asarray(scatter_axis, dtype=float)
    if scatter_axis.shape != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

    n_om = len(ominc)
    neloss = len(eloss)
    gamma_core = _expand_broadening(gamma_c, n_om, 'gamma_c')
    gamma_final = _expand_broadening(gamma_f, neloss, 'gamma_f')

    polarizations = [
        _rixs_polarization_vectors(
            ntrans, thin, thout, phi,
            incoming_kind, alpha, outgoing_kind, beta, scatter_axis,
        )
        for incoming_kind, alpha, outgoing_kind, beta in pol_type
    ]

    # Reusable shifted-system matrix and Krylov solver.
    shifted = hmat_n.duplicate(copy=True)
    ksp = PETSc.KSP().create(hmat_n.getComm())
    ksp.setType(ksp_type)
    ksp.setTolerances(atol=linsys_tol, max_it=linsys_max)
    ksp.setFromOptions()

    rixs = np.zeros((n_om, neloss, len(pol_type)), dtype=float)
    poles = []
    for iom, omega in enumerate(ominc):
        poles_per_om = []
        for ip, (polvec_i, polvec_f) in enumerate(polarizations):
            pole_dict = {'npoles': [], 'eigval': [], 'norm': [],
                         'alpha': [], 'beta': []}
            for ig in range(num_gs):
                z = omega + eval_i[ig] + 1j * gamma_core[iom]

                # A = z*I - H_n
                shifted.zeroEntries()
                shifted.axpy(-1.0, hmat_n,
                             structure=PETSc.Mat.Structure.SAME_NONZERO_PATTERN)
                shifted.shift(z)
                ksp.setOperators(shifted)

                # RHS: D_k |i>  (initial -> intermediate)
                rhs = _apply_trans_combination(trans_op, polvec_i, evec_i[ig])
                sol = rhs.duplicate()
                sol.zeroEntries()
                ksp.solve(rhs, sol)
                if ksp.getConvergedReason() < 0:
                    raise RuntimeError(
                        "KSP did not converge for omega={}, ig={}: reason={}".format(
                            omega, ig, ksp.getConvergedReason()
                        )
                    )

                # Final-state seed: D_k'^dagger x  (intermediate -> initial)
                seed = _apply_trans_combination_adjoint(
                    trans_op, np.conj(polvec_f), sol
                )

                if skip_gs:
                    for gs_vec in evec_i:
                        seed.axpy(-gs_vec.dot(seed), gs_vec)

                alpha_i, beta_i, norm_i = lanczos_tridiagonal(
                    hmat_i, seed, nkryl=nkryl
                )
                pole_dict['npoles'].append(len(alpha_i))
                pole_dict['eigval'].append(eval_i[ig])
                pole_dict['norm'].append(norm_i)
                pole_dict['alpha'].append(alpha_i)
                pole_dict['beta'].append(beta_i)

            poles_per_om.append(pole_dict)
            rixs[iom, :, ip] = get_spectra_from_poles(
                pole_dict, eloss, gamma_final, temperature
            )
        poles.append(poles_per_om)

    return (rixs, poles) if return_poles else rixs
