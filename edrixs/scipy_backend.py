"""SciPy implementation of EDRIXS operator construction and Krylov solvers."""

from __future__ import annotations

from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
import inspect
import multiprocessing as mp
import os
import pickle
import warnings

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, aslinearoperator, gmres, lobpcg

from ._solvers_helpers import _expand_broadening
from .krylov import lanczos_tridiagonal
from .photon_transition import (
    dipole_polvec_xas, dipole_polvec_rixs, quadrupole_polvec, unit_wavevector,
)
from .plot_spectrum import get_spectra_from_poles

__all__ = [
    'owns_operator_scipy', 'owns_operator_dense',
    'get_H_scipy', 'get_transition_operators_scipy',
    'get_H_dense', 'get_transition_operators_dense',
    'ed_scipy', 'xas_scipy', 'rixs_scipy',
    'ed_dense', 'xas_dense', 'rixs_dense',
    'ed_krylov_scipy', 'xas_krylov_scipy',
    'rixs_krylov_scipy', 'rixs_krylov_scipy_parallel',
]


def _backend_kws(backend_kws):
    if backend_kws is None:
        return {}
    if not isinstance(backend_kws, Mapping):
        raise TypeError("backend_kws must be a mapping or None")
    return dict(backend_kws)


def owns_operator_scipy(operator):
    """Return whether ``operator`` is natively accepted by the SciPy backend."""
    return (
        isinstance(operator, np.ndarray)
        or sp.issparse(operator)
        or isinstance(operator, LinearOperator)
    )


def owns_operator_dense(operator):
    """Compatibility recognizer for explicitly requested dense arrays."""
    return isinstance(operator, np.ndarray)


# -----------------------------------------------------------------------------
# SciPy CSR many-body operator construction
# -----------------------------------------------------------------------------


def _count_occupied_before(state, orbital, norbs):
    return (int(state) >> (norbs - orbital)).bit_count()


def two_fermion_csr(emat, left_basis, right_basis=None, tol=1e-10):
    if right_basis is None:
        right_basis = left_basis

    if right_basis.norbs != left_basis.norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    emat = np.asarray(emat)
    norbs = left_basis.norbs
    if emat.shape != (norbs, norbs):
        raise ValueError(
            "emat has shape {}, expected {}".format(
                emat.shape, (norbs, norbs)
            )
        )

    rows, cols, data = [], [], []
    iorb_all, jorb_all = np.nonzero(np.abs(emat) > tol)

    for iorb, jorb in zip(iorb_all, jorb_all):
        bit_j = 1 << (norbs - 1 - jorb)
        bit_i = 1 << (norbs - 1 - iorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_j:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, jorb, norbs)
            state ^= bit_j

            if state & bit_i:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, iorb, norbs)
            state |= bit_i

            try:
                row = left_basis.encode(state)
            except KeyError:
                continue

            rows.append(row)
            cols.append(column)
            data.append(emat[iorb, jorb] * sign_1 * sign_2)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(left_basis), len(right_basis)),
        dtype=np.complex128,
    ).tocsr()


def four_fermion_csr(umat, left_basis, right_basis=None, tol=1e-10):
    if right_basis is None:
        right_basis = left_basis

    norbs = left_basis.norbs
    if right_basis.norbs != norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    umat = np.asarray(umat)
    expected = (norbs, norbs, norbs, norbs)
    if umat.shape != expected:
        raise ValueError(
            "dense umat has shape {}, expected {}".format(umat.shape, expected)
        )

    rows, cols, data = [], [], []
    nonzero = zip(*np.nonzero(np.abs(umat) > tol))

    for lorb, korb, jorb, iorb in nonzero:
        if iorb == jorb or korb == lorb:
            continue

        bit_i = 1 << (norbs - 1 - iorb)
        bit_j = 1 << (norbs - 1 - jorb)
        bit_k = 1 << (norbs - 1 - korb)
        bit_l = 1 << (norbs - 1 - lorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_i:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, iorb, norbs)
            state ^= bit_i

            if not state & bit_j:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, jorb, norbs)
            state ^= bit_j

            if state & bit_k:
                continue
            sign_3 = (-1) ** _count_occupied_before(state, korb, norbs)
            state |= bit_k

            if state & bit_l:
                continue
            sign_4 = (-1) ** _count_occupied_before(state, lorb, norbs)
            state |= bit_l

            try:
                row = left_basis.encode(state)
            except KeyError:
                continue

            rows.append(row)
            cols.append(column)
            data.append(
                umat[lorb, korb, jorb, iorb]
                * sign_1 * sign_2 * sign_3 * sign_4
            )

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(left_basis), len(right_basis)),
        dtype=np.complex128,
    ).tocsr()


def _four_fermion_csr_from_sparse_umat(
        umat, left_basis, right_basis=None, tol=1e-10):
    if right_basis is None:
        right_basis = left_basis

    norbs = left_basis.norbs
    if right_basis.norbs != norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    expected_shape = (norbs * norbs, norbs * norbs)
    if umat.shape != expected_shape:
        raise ValueError(
            "sparse umat has shape {}, expected {}".format(
                umat.shape, expected_shape
            )
        )

    rows, cols, data = [], [], []
    umat = umat.tocoo()

    for flattened_row, flattened_col, value in zip(
            umat.row, umat.col, umat.data):
        if abs(value) <= tol:
            continue

        lorb, korb = divmod(flattened_row, norbs)
        jorb, iorb = divmod(flattened_col, norbs)
        if iorb == jorb or korb == lorb:
            continue

        bit_i = 1 << (norbs - 1 - iorb)
        bit_j = 1 << (norbs - 1 - jorb)
        bit_k = 1 << (norbs - 1 - korb)
        bit_l = 1 << (norbs - 1 - lorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_i:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, iorb, norbs)
            state ^= bit_i

            if not state & bit_j:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, jorb, norbs)
            state ^= bit_j

            if state & bit_k:
                continue
            sign_3 = (-1) ** _count_occupied_before(state, korb, norbs)
            state |= bit_k

            if state & bit_l:
                continue
            sign_4 = (-1) ** _count_occupied_before(state, lorb, norbs)
            state |= bit_l

            try:
                row = left_basis.encode(state)
            except KeyError:
                continue

            rows.append(row)
            cols.append(column)
            data.append(value * sign_1 * sign_2 * sign_3 * sign_4)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(left_basis), len(right_basis)),
        dtype=np.complex128,
    ).tocsr()


def four_fermion_csr_auto(umat, basis, right_basis=None, tol=1e-10):
    if sp.issparse(umat):
        return _four_fermion_csr_from_sparse_umat(
            umat, basis, right_basis=right_basis, tol=tol
        )
    return four_fermion_csr(
        umat, basis, right_basis=right_basis, tol=tol
    )


def get_H_scipy(emat, umat, basis, *, backend_kws=None):
    """Build and return a SciPy CSR many-body Hamiltonian."""
    kws = _backend_kws(backend_kws)
    tol = kws.pop('tol', 1e-10)
    if kws:
        raise TypeError("Unknown SciPy operator-construction options: {}".format(
            sorted(kws)
        ))

    hmat = two_fermion_csr(emat, basis, basis, tol=tol)
    hmat = hmat + four_fermion_csr_auto(umat, basis, tol=tol)
    return hmat.tocsr()


def get_transition_operators_scipy(
        trans_mat, basis_n, basis_i, *, backend_kws=None):
    """Build transition operators as SciPy CSR matrices."""
    kws = _backend_kws(backend_kws)
    tol = kws.pop('tol', 1e-10)
    if kws:
        raise TypeError("Unknown SciPy operator-construction options: {}".format(
            sorted(kws)
        ))

    trans_mat = np.asarray(trans_mat)
    if trans_mat.ndim != 3:
        raise ValueError("trans_mat must be a three-dimensional array")
    return [
        two_fermion_csr(component, basis_n, basis_i, tol=tol)
        for component in trans_mat
    ]


def get_H_dense(emat, umat, basis, *, backend_kws=None):
    """Compatibility dense constructor implemented through SciPy CSR."""
    return get_H_scipy(
        emat, umat, basis, backend_kws=backend_kws
    ).toarray()


def get_transition_operators_dense(
        trans_mat, basis_n, basis_i, *, backend_kws=None):
    """Compatibility dense transition construction through SciPy CSR."""
    return [
        operator.toarray()
        for operator in get_transition_operators_scipy(
            trans_mat, basis_n, basis_i, backend_kws=backend_kws
        )
    ]


# -----------------------------------------------------------------------------
# Backend dispatch entry points
# -----------------------------------------------------------------------------


def ed_scipy(hmat_i, num_evals=1, *, backend_kws=None):
    kws = _backend_kws(backend_kws)
    return ed_krylov_scipy(hmat_i, num_gs=num_evals, **kws)


def xas_scipy(eval_i, evec_i, hmat_n, trans_op, ominc, *,
              gamma_c=0.1, thin=1.0, phi=0.0, pol_type=None,
              temperature=1.0, scatter_axis=None, backend_kws=None):
    kws = _backend_kws(backend_kws)
    return xas_krylov_scipy(
        eval_i, evec_i, hmat_n, trans_op, ominc,
        gamma_c=gamma_c, thin=thin, phi=phi, pol_type=pol_type,
        temperature=temperature, scatter_axis=scatter_axis, **kws
    )


def rixs_scipy(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
               gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
               pol_type=None, temperature=1.0, scatter_axis=None,
               skip_gs=False, return_poles=False, backend_kws=None):
    kws = _backend_kws(backend_kws)
    parallel = bool(kws.pop('parallel', False))
    implementation = (
        rixs_krylov_scipy_parallel if parallel else rixs_krylov_scipy
    )
    return implementation(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
        gamma_c=gamma_c, gamma_f=gamma_f, thin=thin, thout=thout,
        phi=phi, pol_type=pol_type, temperature=temperature,
        scatter_axis=scatter_axis, skip_gs=skip_gs,
        return_poles=return_poles, **kws
    )


# ``dense`` remains an explicit compatibility alias. Dense NumPy arrays are
# still solved by SciPy until a separate NumPy backend is implemented.
ed_dense = ed_scipy
xas_dense = xas_scipy
rixs_dense = rixs_scipy


# -----------------------------------------------------------------------------
# SciPy Krylov implementations
# -----------------------------------------------------------------------------


def ed_krylov_scipy(hmat_i, num_gs=1, blocksize=None, *,
                     tol=1e-10, maxiter=200, seed=None, initial_guess=None,
                     suppress_lobpcg_warnings=True):
    """Compute the lowest retained eigenpairs with SciPy LOBPCG."""
    hmat_i = aslinearoperator(hmat_i)
    if hmat_i.shape[0] != hmat_i.shape[1]:
        raise ValueError("hmat_i must be square")

    num_gs = int(num_gs)
    if num_gs < 1:
        raise ValueError("num_gs must be a positive integer")

    dimension = hmat_i.shape[0]
    if blocksize is None:
        blocksize = num_gs
    blocksize = int(blocksize)
    if blocksize < num_gs:
        raise ValueError("blocksize must be greater than or equal to num_gs")
    if blocksize > dimension:
        raise ValueError("blocksize cannot exceed hmat_i.shape[0]")

    if initial_guess is None:
        rng = np.random.default_rng(seed)
        initial = rng.normal(size=(dimension, blocksize))
    else:
        initial = np.asarray(initial_guess)
        if initial.shape != (dimension, blocksize):
            raise ValueError(
                "initial_guess must have shape {}, got {}".format(
                    (dimension, blocksize), initial.shape
                )
            )

    def solve():
        return lobpcg(
            hmat_i, initial, largest=False, tol=tol, maxiter=maxiter
        )

    if suppress_lobpcg_warnings:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', category=UserWarning,
                message=r'Exited at iteration .*',
            )
            warnings.filterwarnings(
                'ignore', category=UserWarning,
                message=r'Exited postprocessing .*',
            )
            eigenvalues, eigenvectors = solve()
    else:
        eigenvalues, eigenvectors = solve()

    order = np.argsort(eigenvalues)
    return (
        np.asarray(eigenvalues)[order[:num_gs]],
        np.asarray(eigenvectors)[:, order[:num_gs]],
    )


def _apply_linear_combination(operators, coefficients, vector):
    result = None
    for coefficient, operator in zip(coefficients, operators):
        if coefficient == 0:
            continue
        term = coefficient * (operator @ vector)
        result = term if result is None else result + term
    if result is None:
        return np.zeros(operators[0].shape[0], dtype=complex)
    return result


def _check_adjoint_action(operator, name):
    try:
        _ = operator.H @ np.zeros(operator.shape[0], dtype=complex)
    except Exception as exc:
        raise TypeError(
            "{} must provide an adjoint action. For a custom LinearOperator, "
            "define rmatvec or _adjoint.".format(name)
        ) from exc


def _xas_poles_from_start_vectors(eval_i, start_vectors, hmat_n, *, nkryl):
    poles = {'eigval': [], 'npoles': [], 'norm': [], 'alpha': [], 'beta': []}
    effective_nkryl = min(int(nkryl), hmat_n.shape[0])
    for energy, start in zip(eval_i, start_vectors):
        start = np.asarray(start, dtype=complex)
        if np.linalg.norm(start) == 0:
            alpha = np.array([0.0], dtype=float)
            beta = np.array([], dtype=float)
            norm = 0.0
        else:
            alpha, beta, norm = lanczos_tridiagonal(
                hmat_n, start, m=effective_nkryl
            )
        poles['eigval'].append(float(energy))
        poles['npoles'].append(len(alpha))
        poles['norm'].append(float(np.real(norm)))
        poles['alpha'].append(np.asarray(alpha))
        poles['beta'].append(np.asarray(beta))
    return poles


def xas_krylov_scipy(
        eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c=0.1,
        thin=1.0, phi=0.0, pol_type=None, temperature=1.0,
        scatter_axis=None, nkryl=200):
    """Calculate XAS with a SciPy Lanczos continued-fraction solver."""
    eval_i = np.asarray(eval_i, dtype=float)
    evec_i = np.asarray(evec_i, dtype=complex)
    ominc = np.asarray(ominc, dtype=float)
    hmat_n = aslinearoperator(hmat_n)
    trans_op = [aslinearoperator(operator) for operator in trans_op]

    if eval_i.ndim != 1:
        raise ValueError("eval_i must be a one-dimensional array")
    if evec_i.ndim != 2:
        raise ValueError("evec_i must be a two-dimensional array")
    if evec_i.shape[1] != len(eval_i):
        raise ValueError("evec_i columns must correspond to eval_i")
    if hmat_n.shape[0] != hmat_n.shape[1]:
        raise ValueError("hmat_n must be square")

    dim_i, dim_n = evec_i.shape[0], hmat_n.shape[0]
    ntrans = len(trans_op)
    if ntrans not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )
    for index, operator in enumerate(trans_op):
        if operator.shape != (dim_n, dim_i):
            raise ValueError(
                "trans_op[{}] has shape {}, expected {}".format(
                    index, operator.shape, (dim_n, dim_i)
                )
            )

    if pol_type is None:
        pol_type = [('isotropic', 0.0)]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.asarray(scatter_axis, dtype=float)
    if scatter_axis.shape != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

    nkryl = int(nkryl)
    if nkryl < 1:
        raise ValueError("nkryl must be a positive integer")

    gamma_core = _expand_broadening(gamma_c, len(ominc), 'gamma_c')
    spectrum = np.zeros((len(ominc), len(pol_type)), dtype=float)
    wavevector = unit_wavevector(thin, phi, scatter_axis, direction='in')

    for polarization_index, (kind, alpha) in enumerate(pol_type):
        kind = kind.strip().lower()
        if kind == 'isotropic':
            for component in range(ntrans):
                starts = [
                    trans_op[component] @ evec_i[:, state]
                    for state in range(len(eval_i))
                ]
                poles = _xas_poles_from_start_vectors(
                    eval_i, starts, hmat_n, nkryl=nkryl
                )
                spectrum[:, polarization_index] += get_spectra_from_poles(
                    poles, ominc, gamma_core, temperature
                ) / ntrans
            continue

        if kind not in ('linear', 'left', 'right'):
            raise ValueError("Unknown XAS polarization type: {}".format(kind))

        dipole_vector = dipole_polvec_xas(
            thin, phi, alpha, scatter_axis, kind
        )
        polarization_vector = (
            np.asarray(dipole_vector, dtype=complex)
            if ntrans == 3
            else quadrupole_polvec(dipole_vector, wavevector)
        )
        starts = [
            _apply_linear_combination(
                trans_op, polarization_vector, evec_i[:, state]
            )
            for state in range(len(eval_i))
        ]
        poles = _xas_poles_from_start_vectors(
            eval_i, starts, hmat_n, nkryl=nkryl
        )
        spectrum[:, polarization_index] = get_spectra_from_poles(
            poles, ominc, gamma_core, temperature
        )

    return spectrum


def _gmres_scipy_compat(operator, rhs, *, tol, restart, maxiter):
    signature = inspect.signature(gmres)
    if 'rtol' in signature.parameters:
        return gmres(
            operator, rhs, rtol=tol, atol=0.0,
            restart=restart, maxiter=maxiter,
        )
    return gmres(
        operator, rhs, tol=tol, restart=restart, maxiter=maxiter
    )


def _rixs_polarization_vectors(
        ntrans, thin, thout, phi, incoming_kind, alpha,
        outgoing_kind, beta, scatter_axis):
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
        incoming_vector[:] = quadrupole_polvec(
            incoming, incoming_wavevector
        )
        outgoing_vector[:] = quadrupole_polvec(
            outgoing, outgoing_wavevector
        )
    else:
        raise ValueError("ntrans must be 3 or 5")
    return incoming_vector, outgoing_vector


def _rixs_krylov_one_contribution_scipy(
        *, hmat_i, hmat_n, trans_op_H, polvec_f, eval_i, evec_i,
        istate, omega, gamma_c, rhs, skip_gs, nkryl, linsys_tol,
        linsys_maxiter, linsys_restart):
    initial_energy = eval_i[istate]
    if np.linalg.norm(rhs) == 0:
        return {
            'eigval': initial_energy,
            'npoles': 1,
            'norm': 0.0,
            'alpha': np.array([0.0], dtype=float),
            'beta': np.array([], dtype=float),
        }

    shift = omega + initial_energy + 1j * gamma_c
    linear_system = LinearOperator(
        shape=hmat_n.shape,
        matvec=lambda vector, z=shift: z * vector - hmat_n @ vector,
        dtype=np.complex128,
    )
    solution, info = _gmres_scipy_compat(
        linear_system, rhs, tol=linsys_tol,
        restart=linsys_restart, maxiter=linsys_maxiter,
    )
    if info != 0:
        raise RuntimeError(
            "GMRES did not converge for istate={}, omega={}; info={}".format(
                istate, omega, info
            )
        )

    final_vector = _apply_linear_combination(
        trans_op_H, np.conj(polvec_f), solution
    )
    if skip_gs:
        for state in range(len(eval_i)):
            ground_vector = evec_i[:, state]
            final_vector -= ground_vector * np.vdot(
                ground_vector, final_vector
            )

    if np.linalg.norm(final_vector) == 0:
        alpha = np.array([0.0], dtype=float)
        beta = np.array([], dtype=float)
        norm = 0.0
    else:
        alpha, beta, norm = lanczos_tridiagonal(
            hmat_i, final_vector, m=min(int(nkryl), hmat_i.shape[0])
        )

    return {
        'eigval': initial_energy,
        'npoles': len(alpha),
        'norm': norm,
        'alpha': alpha,
        'beta': beta,
    }


def _prepare_rixs(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
        gamma_c, gamma_f, thin, thout, phi, pol_type, scatter_axis,
        skip_gs, nkryl, linsys_tol, linsys_maxiter, linsys_restart):
    eval_i = np.asarray(eval_i, dtype=float)
    evec_i = np.asarray(evec_i, dtype=complex)
    ominc = np.asarray(ominc, dtype=float)
    eloss = np.asarray(eloss, dtype=float)
    hmat_i = aslinearoperator(hmat_i)
    hmat_n = aslinearoperator(hmat_n)
    trans_op = [aslinearoperator(operator) for operator in trans_op]

    if hmat_i.shape[0] != hmat_i.shape[1]:
        raise ValueError("hmat_i must be square")
    if hmat_n.shape[0] != hmat_n.shape[1]:
        raise ValueError("hmat_n must be square")
    if eval_i.ndim != 1:
        raise ValueError("eval_i must be a one-dimensional array")
    if evec_i.ndim != 2:
        raise ValueError("evec_i must be a two-dimensional array")
    if len(eval_i) != evec_i.shape[1]:
        raise ValueError("len(eval_i) must equal evec_i.shape[1]")
    if evec_i.shape[0] != hmat_i.shape[0]:
        raise ValueError("evec_i.shape[0] must equal hmat_i.shape[0]")

    dim_i, dim_n = hmat_i.shape[0], hmat_n.shape[0]
    ntrans = len(trans_op)
    if ntrans not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )
    for index, operator in enumerate(trans_op):
        if operator.shape != (dim_n, dim_i):
            raise ValueError(
                "trans_op[{}] has shape {}, expected {}".format(
                    index, operator.shape, (dim_n, dim_i)
                )
            )
        _check_adjoint_action(operator, 'trans_op[{}]'.format(index))

    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.asarray(scatter_axis, dtype=float)
    if scatter_axis.shape != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

    nkryl = int(nkryl)
    linsys_maxiter = int(linsys_maxiter)
    linsys_restart = int(linsys_restart)
    if nkryl < 1:
        raise ValueError("nkryl must be a positive integer")
    if linsys_maxiter < 1:
        raise ValueError("linsys_maxiter must be a positive integer")
    if linsys_restart < 1:
        raise ValueError("linsys_restart must be a positive integer")

    gamma_core = _expand_broadening(gamma_c, len(ominc), 'gamma_c')
    gamma_final = _expand_broadening(gamma_f, len(eloss), 'gamma_f')
    polarizations = [
        _rixs_polarization_vectors(
            ntrans, thin, thout, phi,
            incoming_kind, alpha, outgoing_kind, beta, scatter_axis,
        )
        for incoming_kind, alpha, outgoing_kind, beta in pol_type
    ]

    return {
        'eval_i': eval_i,
        'evec_i': evec_i,
        'hmat_i': hmat_i,
        'hmat_n': hmat_n,
        'trans_op': trans_op,
        'ominc': ominc,
        'eloss': eloss,
        'gamma_core': gamma_core,
        'gamma_final': gamma_final,
        'polarizations': polarizations,
        'temperature': None,
        'skip_gs': bool(skip_gs),
        'nkryl': nkryl,
        'linsys_tol': float(linsys_tol),
        'linsys_maxiter': linsys_maxiter,
        'linsys_restart': linsys_restart,
    }


def _pole_dict_from_records(records):
    keys = ('npoles', 'eigval', 'norm', 'alpha', 'beta')
    return {key: [record[key] for record in records] for key in keys}


def _compute_rixs_records_serial(problem):
    eval_i = problem['eval_i']
    evec_i = problem['evec_i']
    trans_op = problem['trans_op']
    trans_op_H = [operator.H for operator in trans_op]
    records = [
        [
            [None for _ in range(len(eval_i))]
            for _ in range(len(problem['polarizations']))
        ]
        for _ in range(len(problem['ominc']))
    ]

    for incident_index, omega in enumerate(problem['ominc']):
        for polarization_index, (polvec_i, polvec_f) in enumerate(
                problem['polarizations']):
            for initial_index in range(len(eval_i)):
                rhs = _apply_linear_combination(
                    trans_op, polvec_i, evec_i[:, initial_index]
                )
                records[incident_index][polarization_index][initial_index] = (
                    _rixs_krylov_one_contribution_scipy(
                        hmat_i=problem['hmat_i'],
                        hmat_n=problem['hmat_n'],
                        trans_op_H=trans_op_H,
                        polvec_f=polvec_f,
                        eval_i=eval_i,
                        evec_i=evec_i,
                        istate=initial_index,
                        omega=omega,
                        gamma_c=problem['gamma_core'][incident_index],
                        rhs=rhs,
                        skip_gs=problem['skip_gs'],
                        nkryl=problem['nkryl'],
                        linsys_tol=problem['linsys_tol'],
                        linsys_maxiter=problem['linsys_maxiter'],
                        linsys_restart=problem['linsys_restart'],
                    )
                )
    return records


def _assemble_rixs(problem, records, temperature, return_poles):
    spectrum = np.zeros(
        (
            len(problem['ominc']),
            len(problem['eloss']),
            len(problem['polarizations']),
        ),
        dtype=float,
    )
    poles_all = [
        [None for _ in problem['polarizations']]
        for _ in problem['ominc']
    ]
    for incident_index in range(len(problem['ominc'])):
        for polarization_index in range(len(problem['polarizations'])):
            poles = _pole_dict_from_records(
                records[incident_index][polarization_index]
            )
            poles_all[incident_index][polarization_index] = poles
            spectrum[incident_index, :, polarization_index] = (
                get_spectra_from_poles(
                    poles,
                    problem['eloss'],
                    problem['gamma_final'],
                    temperature,
                )
            )
    return (spectrum, poles_all) if return_poles else spectrum


def rixs_krylov_scipy(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
        gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
        pol_type=None, temperature=1.0, scatter_axis=None,
        skip_gs=False, nkryl=200, linsys_tol=1e-9,
        linsys_maxiter=50000, linsys_restart=200, return_poles=False):
    """Calculate RIXS with Python-native serial loops."""
    problem = _prepare_rixs(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
        gamma_c=gamma_c, gamma_f=gamma_f, thin=thin, thout=thout,
        phi=phi, pol_type=pol_type, scatter_axis=scatter_axis,
        skip_gs=skip_gs, nkryl=nkryl, linsys_tol=linsys_tol,
        linsys_maxiter=linsys_maxiter, linsys_restart=linsys_restart,
    )
    records = _compute_rixs_records_serial(problem)
    return _assemble_rixs(problem, records, temperature, return_poles)


_RIXS_WORKER_STATE = None
_RIXS_THREADPOOL_LIMITER = None


def _available_cpu_count():
    process_cpu_count = getattr(os, 'process_cpu_count', None)
    if process_cpu_count is not None:
        count = process_cpu_count()
        if count is not None:
            return max(1, int(count))
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, int(os.cpu_count() or 1))


def _limit_worker_native_threads(nthreads):
    global _RIXS_THREADPOOL_LIMITER
    for variable in (
        'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
    ):
        os.environ[variable] = str(int(nthreads))
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        _RIXS_THREADPOOL_LIMITER = None
    else:
        _RIXS_THREADPOOL_LIMITER = threadpool_limits(limits=int(nthreads))


def _operator_for_process_pool(operator, name):
    if sp.issparse(operator) or isinstance(operator, np.ndarray):
        candidate = operator
    else:
        candidate = getattr(operator, 'A', operator)
    try:
        pickle.dumps(candidate, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:
        raise TypeError(
            "{} is not serializable for process-based parallel RIXS. "
            "Use a SciPy sparse matrix, a NumPy array, a "
            "MatrixLinearOperator, or use rixs_krylov_scipy.".format(name)
        ) from exc
    return candidate


def _init_rixs_worker(
        hmat_i, hmat_n, trans_op, eval_i, evec_i, ominc, gamma_core,
        polarizations, skip_gs, nkryl, linsys_tol, linsys_maxiter,
        linsys_restart, blas_threads):
    global _RIXS_WORKER_STATE
    _limit_worker_native_threads(blas_threads)
    trans_op = [aslinearoperator(operator) for operator in trans_op]
    _RIXS_WORKER_STATE = {
        'hmat_i': aslinearoperator(hmat_i),
        'hmat_n': aslinearoperator(hmat_n),
        'trans_op': trans_op,
        'trans_op_H': [operator.H for operator in trans_op],
        'eval_i': np.asarray(eval_i, dtype=float),
        'evec_i': np.asarray(evec_i, dtype=complex),
        'ominc': np.asarray(ominc, dtype=float),
        'gamma_core': np.asarray(gamma_core, dtype=float),
        'polarizations': polarizations,
        'skip_gs': bool(skip_gs),
        'nkryl': int(nkryl),
        'linsys_tol': float(linsys_tol),
        'linsys_maxiter': int(linsys_maxiter),
        'linsys_restart': int(linsys_restart),
    }


def _rixs_pool_job(job):
    incident_index, polarization_index, initial_index = job
    state = _RIXS_WORKER_STATE
    if state is None:
        raise RuntimeError("RIXS worker was not initialized")

    polvec_i, polvec_f = state['polarizations'][polarization_index]
    rhs = _apply_linear_combination(
        state['trans_op'], polvec_i, state['evec_i'][:, initial_index]
    )
    record = _rixs_krylov_one_contribution_scipy(
        hmat_i=state['hmat_i'], hmat_n=state['hmat_n'],
        trans_op_H=state['trans_op_H'], polvec_f=polvec_f,
        eval_i=state['eval_i'], evec_i=state['evec_i'],
        istate=initial_index, omega=state['ominc'][incident_index],
        gamma_c=state['gamma_core'][incident_index], rhs=rhs,
        skip_gs=state['skip_gs'], nkryl=state['nkryl'],
        linsys_tol=state['linsys_tol'],
        linsys_maxiter=state['linsys_maxiter'],
        linsys_restart=state['linsys_restart'],
    )
    return incident_index, polarization_index, initial_index, record


def rixs_krylov_scipy_parallel(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
        gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
        pol_type=None, temperature=1.0, scatter_axis=None,
        skip_gs=False, nkryl=200, linsys_tol=1e-9,
        linsys_maxiter=50000, linsys_restart=200, return_poles=False,
        workers=None, blas_threads=1, mp_start_method=None):
    """Calculate RIXS with a process pool over energy/polarization/state jobs."""
    problem = _prepare_rixs(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
        gamma_c=gamma_c, gamma_f=gamma_f, thin=thin, thout=thout,
        phi=phi, pol_type=pol_type, scatter_axis=scatter_axis,
        skip_gs=skip_gs, nkryl=nkryl, linsys_tol=linsys_tol,
        linsys_maxiter=linsys_maxiter, linsys_restart=linsys_restart,
    )

    blas_threads = int(blas_threads)
    if blas_threads < 1:
        raise ValueError("blas_threads must be a positive integer")

    jobs = [
        (incident_index, polarization_index, initial_index)
        for incident_index in range(len(problem['ominc']))
        for polarization_index in range(len(problem['polarizations']))
        for initial_index in range(len(problem['eval_i']))
    ]
    if workers is None:
        workers = _available_cpu_count()
    workers = int(workers)
    if workers < 1:
        raise ValueError("workers must be None or a positive integer")
    workers = min(workers, max(1, len(jobs)))

    if workers == 1 or len(jobs) <= 1:
        records = _compute_rixs_records_serial(problem)
        return _assemble_rixs(problem, records, temperature, return_poles)

    pool_hmat_i = _operator_for_process_pool(problem['hmat_i'], 'hmat_i')
    pool_hmat_n = _operator_for_process_pool(problem['hmat_n'], 'hmat_n')
    pool_trans_op = [
        _operator_for_process_pool(operator, 'trans_op[{}]'.format(index))
        for index, operator in enumerate(problem['trans_op'])
    ]
    context = (
        mp.get_context() if mp_start_method is None
        else mp.get_context(mp_start_method)
    )
    records = [
        [
            [None for _ in problem['eval_i']]
            for _ in problem['polarizations']
        ]
        for _ in problem['ominc']
    ]

    with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=_init_rixs_worker,
            initargs=(
                pool_hmat_i, pool_hmat_n, pool_trans_op,
                problem['eval_i'], problem['evec_i'], problem['ominc'],
                problem['gamma_core'], problem['polarizations'],
                problem['skip_gs'], problem['nkryl'],
                problem['linsys_tol'], problem['linsys_maxiter'],
                problem['linsys_restart'], blas_threads,
            )) as executor:
        future_to_job = {
            executor.submit(_rixs_pool_job, job): job for job in jobs
        }
        for future in as_completed(future_to_job):
            incident_index, polarization_index, initial_index = future_to_job[future]
            try:
                _, _, _, record = future.result()
            except Exception as exc:
                for pending in future_to_job:
                    pending.cancel()
                raise RuntimeError(
                    "Parallel RIXS job failed for incident index {}, "
                    "polarization index {}, initial-state index {}".format(
                        incident_index, polarization_index, initial_index
                    )
                ) from exc
            records[incident_index][polarization_index][initial_index] = record

    return _assemble_rixs(problem, records, temperature, return_poles)
