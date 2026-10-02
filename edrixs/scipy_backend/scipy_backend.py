"""SciPy implementation of EDRIXS operator construction and Krylov solvers."""

from __future__ import annotations

import inspect
import warnings

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, aslinearoperator, gmres, lobpcg

from .._solvers_helpers import _expand_broadening
from .krylov import lanczos_tridiagonal
from .options import OPTIONS, validate_options
from ..photon_transition import (
    dipole_polvec_xas, dipole_polvec_rixs, quadrupole_polvec, unit_wavevector,
)
from ..poles import get_spectra_from_poles

__all__ = [
    'owns_operator_scipy',
    'build_op_scipy',
    'ed_scipy', 'xas_scipy', 'rixs_scipy',
]


def ed_scipy(hmat_i, num_evals=1, *, shift=0.0, backend_kws=None):
    """
    Compute the lowest initial-state eigenpairs using SciPy LOBPCG.

    Parameters
    ----------
    hmat_i : ndarray, sparse matrix or LinearOperator
        Square Hermitian initial/final Hamiltonian.
    num_evals : int, optional
        Number of lowest-energy eigenpairs to return.
    shift : float, optional
        Real energy offset subtracted before diagonalization and restored to
        the returned eigenvalues. The input Hamiltonian is not modified.
    backend_kws : mapping, optional
        LOBPCG options: ``blocksize`` (defaults to ``num_evals``), ``tol``,
        ``maxiter``, ``seed``, ``initial_guess``, and
        ``suppress_lobpcg_warnings``. An initial guess must have shape
        ``(dim_i, blocksize)``. Extra eigenpairs from a larger block are
        discarded after sorting.

    Returns
    -------
    eval_i : ndarray
        Lowest eigenvalues on the original energy reference, shape
        ``(num_evals,)``.
    evec_i : ndarray
        Corresponding eigenvectors, shape ``(dim_i, num_evals)``.
    """
    kws = {name: option.default for name, option in OPTIONS['ed'].items()}
    kws.update(validate_options('ed', backend_kws))
    blocksize = kws['blocksize']
    if shift != 0:
        hmat_i = _shift_hamiltonian(hmat_i, shift)
    hmat_i = aslinearoperator(hmat_i)
    if hmat_i.shape[0] != hmat_i.shape[1]:
        raise ValueError("hmat_i must be square")

    if num_evals < 1:
        raise ValueError("num_evals must be a positive integer")

    dimension = hmat_i.shape[0]
    if blocksize is None:
        blocksize = num_evals
    if blocksize < num_evals:
        raise ValueError("blocksize must be greater than or equal to num_evals")
    if blocksize > dimension:
        raise ValueError("blocksize cannot exceed hmat_i.shape[0]")

    initial_guess = kws['initial_guess']
    if initial_guess is None:
        rng = np.random.default_rng(kws['seed'])
        initial = rng.normal(size=(dimension, blocksize))
    else:
        initial = np.asarray(initial_guess)
        if initial.shape != (dimension, blocksize):
            raise ValueError(
                "initial_guess must have shape {}, got {}".format(
                    (dimension, blocksize), initial.shape
                )
            )

    with warnings.catch_warnings():
        if kws['suppress_lobpcg_warnings']:
            warnings.filterwarnings(
                'ignore', category=UserWarning,
                message=r'Exited at iteration .*',
            )
            warnings.filterwarnings(
                'ignore', category=UserWarning,
                message=r'Exited postprocessing .*',
            )
        eigenvalues, eigenvectors = lobpcg(
            hmat_i, initial, largest=False, tol=kws['tol'], maxiter=kws['maxiter']
        )

    order = np.argsort(eigenvalues)[:num_evals]
    return eigenvalues[order] + shift, eigenvectors[:, order]


def xas_scipy(
        eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c=0.1,
        thin=1.0, phi=0.0, pol_type=None, temperature=1.0,
        scatter_axis=None, backend_kws=None):
    """
    Calculate XAS spectra with a SciPy Lanczos continued-fraction solver.

    ``eval_i`` and ``evec_i`` are the retained initial states, normally
    returned by :func:`ed_scipy`.  For every retained initial state the
    transition operator is applied in the original Fock basis, and a Lanczos
    tridiagonalization of the intermediate Hamiltonian generates the pole
    representation consumed by :func:`get_spectra_from_poles`.

    The first incident energy and lowest retained energy set the internal
    reference, keeping the intermediate Hamiltonian near zero on resonance.

    Parameters
    ----------
    eval_i : 1d array
        Energies of the retained initial states.

    evec_i : 2d array
        Retained initial-state eigenvectors. Column ``i`` corresponds to
        ``eval_i[i]``.

    hmat_n : sparse matrix or scipy.sparse.linalg.LinearOperator
        Intermediate-state Hamiltonian.

    trans_op : sequence
        Transition operators mapping the initial Hilbert space to the
        intermediate Hilbert space. The sequence length must be 3 for a
        dipole transition or 5 for a quadrupole transition.

    ominc : 1d array
        Incident photon-energy grid.

    gamma_c : float or 1d array, optional
        Core-hole lifetime broadening. It can be scalar or have the same shape
        as ``ominc``.

    thin, phi : float, optional
        Incoming angle and azimuthal angle, in radians.

    pol_type : sequence of tuple, optional
        Incoming polarizations. Each entry is ``(kind, alpha)`` where ``kind``
        is ``'linear'``, ``'left'``, ``'right'``, or ``'isotropic'``. The
        default is isotropic polarization.

    temperature : float, optional
        Temperature in kelvin used for the Boltzmann weights of ``eval_i``.

    scatter_axis : (3, 3) array, optional
        Scattering coordinate axes. The default is the identity matrix.

    backend_kws : mapping, optional
        ``nkryl`` sets the maximum intermediate-state Lanczos dimension.

    Returns
    -------
    xas : 2d ndarray
        Spectrum with shape ``(len(ominc), len(pol_type))``.
    """
    kws = validate_options('xas', backend_kws)
    nkryl = kws.get('nkryl', OPTIONS['xas']['nkryl'].default)
    eval_i, evec_i, trans_op = _prepare_transition_inputs(
        eval_i, evec_i, trans_op
    )
    ominc = np.asarray(ominc)
    energy_ref = min(eval_i, default=0.0)
    omega_ref = ominc[0] if len(ominc) else 0.0
    hmat_n = aslinearoperator(_shift_hamiltonian(hmat_n, energy_ref, omega_ref))
    eval_i = eval_i - energy_ref
    ominc = ominc - omega_ref
    ntrans = len(trans_op)

    if pol_type is None:
        pol_type = [('isotropic', 0.0)]
    if scatter_axis is not None and np.shape(scatter_axis) != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

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
            dipole_vector
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


def rixs_scipy(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
        gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
        pol_type=None, temperature=1.0, scatter_axis=None,
        skip_gs=False, return_poles=False, backend_kws=None):
    """
    Calculate RIXS spectra with the SciPy Krylov correction-vector solver.

    For each incident energy and polarization, collect the retained-state
    poles and evaluate their thermally weighted spectrum.

    Hamiltonians are centered using the lowest retained energy and the first
    incident energy. Returned poles retain the original absolute energies.

    Parameters
    ----------
    eval_i : 1d array
        Energies of the retained initial states.
    evec_i : 2d array
        Retained initial-state eigenvectors in the basis of ``hmat_i``.
        Column ``i`` corresponds to ``eval_i[i]``.
    hmat_i, hmat_n : sparse matrix or LinearOperator
        Initial/final and intermediate Hamiltonians.
    trans_op : sequence
        Transition operators mapping the initial/final Hilbert space to the
        intermediate Hilbert space. The sequence length must be 3 or 5.
        Custom LinearOperators must supply an adjoint action (``rmatvec``
        or ``_adjoint``) for emission.
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
        If true, omit transitions into the retained initial-state subspace
        from the final-state spectrum.
    return_poles : bool, optional
        Return the nested pole dictionaries together with the spectrum.
    backend_kws : mapping, optional
        ``nkryl`` sets the final-state Lanczos dimension. ``linsys_tol``,
        ``linsys_maxiter``, and ``linsys_restart`` control the intermediate
        GMRES correction-vector solve. ``linsys_tol`` sets the relative
        tolerance (``rtol``, or ``tol`` in older SciPy). When the ``rtol``
        API is available, the absolute tolerance ``atol`` is set to zero.

    Returns
    -------
    rixs : 3d ndarray
        Spectrum with shape ``(len(ominc), len(eloss), len(pol_type))``.
    poles : list, optional
        Returned only when ``return_poles`` is true.
    """
    kws = {name: option.default for name, option in OPTIONS['rixs'].items()}
    kws.update(validate_options('rixs', backend_kws))
    eval_i, evec_i, trans_op = _prepare_transition_inputs(
        eval_i, evec_i, trans_op
    )
    energy_ref = min(eval_i, default=0.0)
    omega_ref = ominc[0] if len(ominc) else 0.0
    hmat_i = aslinearoperator(_shift_hamiltonian(hmat_i, energy_ref))
    hmat_n = aslinearoperator(_shift_hamiltonian(hmat_n, energy_ref, omega_ref))
    eval_relative = eval_i - energy_ref
    trans_op_H = [operator.H for operator in trans_op]
    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if scatter_axis is not None and np.shape(scatter_axis) != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")

    gamma_core = _expand_broadening(gamma_c, len(ominc), 'gamma_c')
    gamma_final = _expand_broadening(gamma_f, len(eloss), 'gamma_f')
    polarizations = []
    if len(trans_op) == 5:
        wavevector_i = unit_wavevector(thin, phi, scatter_axis, direction='in')
        wavevector_f = unit_wavevector(
            thout, phi, scatter_axis, direction='out'
        )
    for incoming_kind, alpha, outgoing_kind, beta in pol_type:
        incoming, outgoing = dipole_polvec_rixs(
            thin, thout, phi, alpha, beta, scatter_axis,
            (incoming_kind, outgoing_kind),
        )
        if len(trans_op) == 5:
            incoming = quadrupole_polvec(incoming, wavevector_i)
            outgoing = quadrupole_polvec(outgoing, wavevector_f)
        polarizations.append((incoming, outgoing))

    gmres_kws = {
        'restart': kws['linsys_restart'], 'maxiter': kws['linsys_maxiter'],
    }
    if 'rtol' in inspect.signature(gmres).parameters:
        gmres_kws.update(rtol=kws['linsys_tol'], atol=0.0)
    else:
        gmres_kws['tol'] = kws['linsys_tol']

    spectrum = np.zeros((len(ominc), len(eloss), len(pol_type)))
    poles_all = []
    for incident_index, omega in enumerate(ominc):
        incident_poles = []
        for polarization_index, (polvec_i, polvec_f) in enumerate(
                polarizations):
            poles = {
                'eigval': [], 'npoles': [], 'norm': [], 'alpha': [], 'beta': [],
            }
            for initial_index, energy in enumerate(eval_relative):
                rhs = _apply_linear_combination(
                    trans_op, polvec_i, evec_i[:, initial_index]
                )
                alpha, beta, norm = np.array([0.0]), np.array([]), 0.0
                if np.linalg.norm(rhs) != 0:
                    shift = (omega - omega_ref) + energy + 1j * gamma_core[incident_index]
                    linear_system = LinearOperator(
                        hmat_n.shape,
                        matvec=lambda v, z=shift: z * v - hmat_n @ v,
                        dtype=complex,
                    )
                    solution, info = gmres(linear_system, rhs, **gmres_kws)
                    if info != 0:
                        raise RuntimeError(
                            "GMRES did not converge for istate={}, omega={}; "
                            "info={}".format(initial_index, omega, info)
                        )

                    final_vector = _apply_linear_combination(
                        trans_op_H, np.conj(polvec_f), solution
                    )
                    if skip_gs:
                        final_vector = final_vector - evec_i @ (
                            evec_i.conj().T @ final_vector
                        )
                    if np.linalg.norm(final_vector) != 0:
                        alpha, beta, norm = lanczos_tridiagonal(
                            hmat_i, final_vector, m=kws['nkryl']
                        )

                poles['eigval'].append(energy)
                poles['npoles'].append(len(alpha))
                poles['norm'].append(norm)
                poles['alpha'].append(alpha)
                poles['beta'].append(beta)

            spectrum[incident_index, :, polarization_index] = (
                get_spectra_from_poles(poles, eloss, gamma_final, temperature)
            )
            if return_poles:
                # Restore absolute energies only after evaluating the spectrum.
                poles['eigval'] = list(eval_i)
                poles['alpha'] = [
                    alpha + energy_ref if norm != 0 else alpha
                    for alpha, norm in zip(poles['alpha'], poles['norm'])
                ]
                incident_poles.append(poles)
        if return_poles:
            poles_all.append(incident_poles)
    return (spectrum, poles_all) if return_poles else spectrum


def owns_operator_scipy(operator):
    """Return whether ``operator`` is natively accepted by the SciPy backend."""
    return (
        isinstance(operator, np.ndarray)
        or sp.issparse(operator)
        or isinstance(operator, LinearOperator)
    )


# -----------------------------------------------------------------------------
# SciPy CSR many-body operator construction
# -----------------------------------------------------------------------------


def _count_occupied_before(state, orbital_bit):
    """
    Count occupied orbitals below the single-bit mask ``orbital_bit``.
    """
    return (state & (orbital_bit - 1)).bit_count()


def two_fermion_csr(emat, left_basis, right_basis=None, tol=1e-10):
    """
    Build a sparse many-body representation of a one-body operator.

    Parameters
    ----------
    emat : 2d array
        Orbital-space coefficients of :math:`c_i^\\dagger c_j`.
    left_basis, right_basis : FockBasis
        Output and input many-electron bases. ``right_basis`` defaults to
        ``left_basis``.
    tol : float, optional
        Ignore coefficients with magnitude not exceeding this threshold.
    """
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
        bit_j = 1 << int(jorb)
        bit_i = 1 << int(iorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_j:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, bit_j)
            state ^= bit_j

            if state & bit_i:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, bit_i)
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
    """
    Build a sparse many-body representation of a dense Coulomb tensor.
    """
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

        bit_i = 1 << int(iorb)
        bit_j = 1 << int(jorb)
        bit_k = 1 << int(korb)
        bit_l = 1 << int(lorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_i:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, bit_i)
            state ^= bit_i

            if not state & bit_j:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, bit_j)
            state ^= bit_j

            if state & bit_k:
                continue
            sign_3 = (-1) ** _count_occupied_before(state, bit_k)
            state |= bit_k

            if state & bit_l:
                continue
            sign_4 = (-1) ** _count_occupied_before(state, bit_l)
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


def four_fermion_csr_auto(umat, basis, right_basis=None, tol=1e-10):
    """Build a sparse two-body operator from dense or sparse Coulomb data."""
    if not sp.issparse(umat):
        return four_fermion_csr(umat, basis, right_basis=right_basis, tol=tol)

    if right_basis is None:
        right_basis = basis

    norbs = basis.norbs
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

        bit_i = 1 << int(iorb)
        bit_j = 1 << int(jorb)
        bit_k = 1 << int(korb)
        bit_l = 1 << int(lorb)

        for column in range(len(right_basis)):
            state = right_basis.decode(column)
            if not state & bit_i:
                continue
            sign_1 = (-1) ** _count_occupied_before(state, bit_i)
            state ^= bit_i

            if not state & bit_j:
                continue
            sign_2 = (-1) ** _count_occupied_before(state, bit_j)
            state ^= bit_j

            if state & bit_k:
                continue
            sign_3 = (-1) ** _count_occupied_before(state, bit_k)
            state |= bit_k

            if state & bit_l:
                continue
            sign_4 = (-1) ** _count_occupied_before(state, bit_l)
            state |= bit_l

            try:
                row = basis.encode(state)
            except KeyError:
                continue

            rows.append(row)
            cols.append(column)
            data.append(value * sign_1 * sign_2 * sign_3 * sign_4)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(basis), len(right_basis)),
        dtype=np.complex128,
    ).tocsr()


def build_op_scipy(emat, umat, lb, rb=None, *, use_numba=False, backend_kws=None):
    """Build and return a SciPy CSR many-body operator."""
    from .hash_basis_methods import build_op_scipy_matrix

    kws = validate_options('build_op', backend_kws)
    tol = kws.get('tol', OPTIONS['build_op']['tol'].default)

    return build_op_scipy_matrix(
        emat, umat, lb, rb, tol=tol, use_numba=use_numba
    )


# -----------------------------------------------------------------------------
# Solver helpers
# -----------------------------------------------------------------------------


def _shift_hamiltonian(hamiltonian, *shifts):
    """Subtract real offsets before matvecs, without modifying the input.

    Subtract offsets separately to preserve a small photon-energy correction
    beside a large common energy. Opaque LinearOperators must shift lazily;
    their original matvec can still lose precision through cancellation.
    """
    if isinstance(hamiltonian, LinearOperator):
        identity = LinearOperator(
            hamiltonian.shape, dtype=float,
            matvec=lambda x: x, rmatvec=lambda x: x,
            matmat=lambda x: x, rmatmat=lambda x: x,
        )
        for shift in shifts:
            hamiltonian = hamiltonian - shift * identity
        return hamiltonian
    if sp.issparse(hamiltonian):
        result = sp.csr_matrix(
            hamiltonian, dtype=np.result_type(hamiltonian.dtype, float), copy=True
        )
        diagonal = result.diagonal()
        for shift in shifts:
            diagonal -= shift
        result.setdiag(diagonal)
        return result
    hamiltonian = np.asarray(hamiltonian)
    result = np.array(
        hamiltonian, dtype=np.result_type(hamiltonian.dtype, float), copy=True
    )
    indices = np.diag_indices_from(result)
    for shift in shifts:
        result[indices] -= shift
    return result


def _apply_linear_combination(operators, coefficients, vector):
    """
    Apply sum_i coeffs[i] ops[i] to vec without constructing the summed operator.
    """
    result = None
    for coefficient, operator in zip(coefficients, operators):
        if coefficient == 0:
            continue
        term = coefficient * (operator @ vector)
        result = term if result is None else result + term
    if result is None:
        return np.zeros(operators[0].shape[0], dtype=complex)
    return result


def _prepare_transition_inputs(eval_i, evec_i, trans_op):
    """Normalize XAS/RIXS inputs and check state and component counts."""
    eval_i = np.asarray(eval_i)
    evec_i = np.asarray(evec_i)
    trans_op = [aslinearoperator(operator) for operator in trans_op]

    if eval_i.ndim != 1:
        raise ValueError("eval_i must be a one-dimensional array")
    if evec_i.ndim != 2:
        raise ValueError("evec_i must be a two-dimensional array")
    if evec_i.shape[1] != len(eval_i):
        raise ValueError("evec_i columns must correspond to eval_i")
    if len(trans_op) not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )
    return eval_i, evec_i, trans_op


def _xas_poles_from_start_vectors(eval_i, start_vectors, hmat_n, *, nkryl):
    """
    Build an EDRIXS-compatible XAS pole dictionary.
    """
    poles = {'eigval': [], 'npoles': [], 'norm': [], 'alpha': [], 'beta': []}
    for energy, start in zip(eval_i, start_vectors):
        if np.linalg.norm(start) == 0:
            alpha = np.array([0.0], dtype=float)
            beta = np.array([], dtype=float)
            norm = 0.0
        else:
            alpha, beta, norm = lanczos_tridiagonal(
                hmat_n, start, m=nkryl
            )
        poles['eigval'].append(energy)
        poles['npoles'].append(len(alpha))
        poles['norm'].append(norm)
        poles['alpha'].append(alpha)
        poles['beta'].append(beta)
    return poles
