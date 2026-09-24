"""Exact dense-matrix implementation of the staged solver interface."""

from __future__ import annotations

import numpy as np
import scipy.linalg

from .._solvers_helpers import _expand_broadening
from ..photon_transition import (
    dipole_polvec_rixs,
    dipole_polvec_xas,
    quadrupole_polvec,
    unit_wavevector,
)
from ..scipy_backend.krylov import lanczos_tridiagonal
from ..scipy_backend.scipy_backend import build_op_scipy
from ..utils import boltz_dist
from .options import validate_options

__all__ = [
    'owns_operator_dense',
    'build_op_dense',
    'ed_dense',
    'xas_dense',
    'rixs_dense',
]


def owns_operator_dense(operator):
    """Return whether ``operator`` is a dense NumPy matrix."""
    return isinstance(operator, np.ndarray) and operator.ndim == 2


def _dense_matrix(operator, name, *, square=False):
    """Return one validated dense matrix."""
    matrix = np.asarray(operator)
    if matrix.ndim != 2:
        raise ValueError("{} must be a two-dimensional array".format(name))
    if square and matrix.shape[0] != matrix.shape[1]:
        raise ValueError("{} must be square".format(name))
    return matrix


def _transition_array(trans_op, dim_n, dim_i):
    """Return validated dense transition operators."""
    transitions = np.asarray(trans_op, dtype=complex)
    if transitions.ndim != 3:
        raise ValueError("trans_op must be a three-dimensional array")
    if transitions.shape[0] not in (3, 5):
        raise ValueError(
            "len(trans_op) must be 3 for dipole or 5 for quadrupole transitions"
        )
    expected = (transitions.shape[0], dim_n, dim_i)
    if transitions.shape != expected:
        raise ValueError(
            "trans_op has shape {}, expected {}".format(
                transitions.shape, expected
            )
        )
    return transitions


def _initial_states(eval_i, evec_i, dim_i=None):
    """Return validated retained initial eigenpairs."""
    energies = np.asarray(eval_i, dtype=float)
    vectors = np.asarray(evec_i, dtype=complex)
    if energies.ndim != 1:
        raise ValueError("eval_i must be a one-dimensional array")
    if len(energies) == 0:
        raise ValueError("eval_i must contain at least one state")
    if vectors.ndim != 2:
        raise ValueError("evec_i must be a two-dimensional array")
    if dim_i is None:
        dim_i = vectors.shape[0]
    if vectors.shape != (dim_i, len(energies)):
        raise ValueError(
            "evec_i has shape {}, expected {}".format(
                vectors.shape, (dim_i, len(energies))
            )
        )
    return energies, vectors


def _scattering_axis(scatter_axis):
    """Return a validated scattering frame."""
    if scatter_axis is None:
        return np.eye(3)
    axis = np.asarray(scatter_axis, dtype=float)
    if axis.shape != (3, 3):
        raise ValueError("scatter_axis must have shape (3, 3)")
    return axis


def _energy_mesh(values, name):
    """Return a validated one-dimensional energy mesh."""
    mesh = np.asarray(values, dtype=float)
    if mesh.ndim != 1:
        raise ValueError("{} must be a one-dimensional array".format(name))
    return mesh


def build_op_dense(
        emat, umat, lb, rb=None, *, use_numba=False, backend_kws=None):
    """Build a many-body operator and return it as a dense NumPy matrix."""
    kws = validate_options('build_op', backend_kws)
    return build_op_scipy(
        emat, umat, lb, rb, use_numba=use_numba, backend_kws=kws
    ).toarray()


def ed_dense(hmat_i, num_evals=1, *, backend_kws=None):
    """Exactly diagonalize a dense Hermitian Hamiltonian."""
    validate_options('ed', backend_kws)
    hamiltonian = _dense_matrix(hmat_i, 'hmat_i', square=True)
    num_evals = int(num_evals)
    dimension = hamiltonian.shape[0]
    if num_evals < 1:
        raise ValueError("num_evals must be a positive integer")
    if num_evals > dimension:
        raise ValueError("num_evals cannot exceed hmat_i.shape[0]")

    subset = None if num_evals == dimension else (0, num_evals - 1)
    return scipy.linalg.eigh(hamiltonian, subset_by_index=subset)


def _xas_polarization_vector(
        ntrans, kind, alpha, thin, phi, scatter_axis, wavevector):
    """Return one XAS polarization vector in transition-component space."""
    kind = kind.strip().lower()
    if kind not in ('linear', 'left', 'right'):
        raise ValueError("Unknown XAS polarization type: {}".format(kind))
    dipole = dipole_polvec_xas(
        thin, phi, alpha, scatter_axis, kind
    )
    return (
        np.asarray(dipole, dtype=complex)
        if ntrans == 3
        else quadrupole_polvec(dipole, wavevector)
    )


def xas_dense(
        eval_i, evec_i, hmat_n, trans_op, ominc, *, gamma_c=0.1,
        thin=1.0, phi=0.0, pol_type=None, temperature=1.0,
        scatter_axis=None, backend_kws=None):
    """Calculate XAS by summing over the exact intermediate eigenstates."""
    validate_options('xas', backend_kws)
    hamiltonian_n = _dense_matrix(hmat_n, 'hmat_n', square=True)
    energies_i, vectors_i = _initial_states(eval_i, evec_i)
    transitions = _transition_array(
        trans_op, hamiltonian_n.shape[0], vectors_i.shape[0]
    )
    incident_energies = _energy_mesh(ominc, 'ominc')
    axis = _scattering_axis(scatter_axis)
    gamma_core = _expand_broadening(
        gamma_c, len(incident_energies), 'gamma_c'
    )
    if pol_type is None:
        pol_type = [('isotropic', 0.0)]

    energies_n, vectors_n = scipy.linalg.eigh(hamiltonian_n)
    transitions_eigen = np.einsum(
        'an,kab,bi->kni',
        vectors_n.conj(),
        transitions,
        vectors_i,
        optimize=True,
    )
    excitation_energies = energies_n[:, None] - energies_i[None, :]
    lorentzian = (
        gamma_core[:, None, None]
        / np.pi
        / (
            (incident_energies[:, None, None]
             - excitation_energies[None, :, :]) ** 2
            + gamma_core[:, None, None] ** 2
        )
    )
    probabilities = boltz_dist(energies_i, temperature)
    spectrum = np.zeros((len(incident_energies), len(pol_type)), dtype=float)
    wavevector = unit_wavevector(thin, phi, axis, direction='in')
    ntrans = transitions.shape[0]

    for channel, (kind, alpha) in enumerate(pol_type):
        kind = kind.strip().lower()
        if kind == 'isotropic':
            strengths = np.sum(np.abs(transitions_eigen) ** 2, axis=0) / ntrans
        else:
            polarization = _xas_polarization_vector(
                ntrans, kind, alpha, thin, phi, axis, wavevector
            )
            amplitudes = np.einsum(
                'k,kni->ni', polarization, transitions_eigen,
                optimize=True,
            )
            strengths = np.abs(amplitudes) ** 2
        spectrum[:, channel] = np.sum(
            lorentzian * strengths[None, :, :]
            * probabilities[None, None, :],
            axis=(1, 2),
        )
    return spectrum


def _rixs_polarization_vectors(
        ntrans, thin, thout, phi, incoming_kind, alpha,
        outgoing_kind, beta, scatter_axis):
    """Return incoming and outgoing RIXS transition-component vectors."""
    incoming_kind = incoming_kind.strip().lower()
    outgoing_kind = outgoing_kind.strip().lower()
    incoming, outgoing = dipole_polvec_rixs(
        thin, thout, phi, alpha, beta, scatter_axis,
        (incoming_kind, outgoing_kind),
    )
    if ntrans == 3:
        return (
            np.asarray(incoming, dtype=complex),
            np.asarray(outgoing, dtype=complex),
        )
    incoming_wavevector = unit_wavevector(
        thin, phi, scatter_axis, direction='in'
    )
    outgoing_wavevector = unit_wavevector(
        thout, phi, scatter_axis, direction='out'
    )
    return (
        quadrupole_polvec(incoming, incoming_wavevector),
        quadrupole_polvec(outgoing, outgoing_wavevector),
    )


def _pole_record(hamiltonian, vector, initial_energy):
    """Represent one exact dense final-state vector in pole-dictionary form."""
    if np.linalg.norm(vector) == 0:
        alpha = np.array([0.0], dtype=float)
        beta = np.array([], dtype=float)
        norm = 0.0
    else:
        alpha, beta, norm = lanczos_tridiagonal(
            hamiltonian, vector, m=hamiltonian.shape[0]
        )
    return {
        'eigval': float(initial_energy),
        'npoles': len(alpha),
        'norm': float(np.real(norm)),
        'alpha': np.asarray(alpha),
        'beta': np.asarray(beta),
    }


def _pole_dict(records):
    """Combine per-initial-state pole records into the public dictionary."""
    keys = ('npoles', 'eigval', 'norm', 'alpha', 'beta')
    return {key: [record[key] for record in records] for key in keys}


def rixs_dense(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
        gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
        pol_type=None, temperature=1.0, scatter_axis=None,
        skip_gs=False, return_poles=False, backend_kws=None):
    """Calculate RIXS by summing over exact intermediate and final states."""
    validate_options('rixs', backend_kws)
    hamiltonian_i = _dense_matrix(hmat_i, 'hmat_i', square=True)
    hamiltonian_n = _dense_matrix(hmat_n, 'hmat_n', square=True)
    energies_i, vectors_i = _initial_states(
        eval_i, evec_i, hamiltonian_i.shape[0]
    )
    transitions = _transition_array(
        trans_op, hamiltonian_n.shape[0], hamiltonian_i.shape[0]
    )
    incident_energies = _energy_mesh(ominc, 'ominc')
    energy_losses = _energy_mesh(eloss, 'eloss')
    gamma_core = _expand_broadening(
        gamma_c, len(incident_energies), 'gamma_c'
    )
    gamma_final = _expand_broadening(
        gamma_f, len(energy_losses), 'gamma_f'
    )
    axis = _scattering_axis(scatter_axis)
    if pol_type is None:
        pol_type = [('linear', 0.0, 'linear', 0.0)]

    energies_f, vectors_f = scipy.linalg.eigh(hamiltonian_i)
    energies_n, vectors_n = scipy.linalg.eigh(hamiltonian_n)
    absorption = np.einsum(
        'an,kab,bi->kni',
        vectors_n.conj(), transitions, vectors_i, optimize=True,
    )
    emission = np.einsum(
        'af,kab,bn->kfn',
        vectors_f.conj(), transitions.conj().transpose(0, 2, 1),
        vectors_n, optimize=True,
    )
    final_energy_changes = energies_f[:, None] - energies_i[None, :]
    final_lorentzian = (
        gamma_final[:, None, None]
        / np.pi
        / (
            (energy_losses[:, None, None]
             - final_energy_changes[None, :, :]) ** 2
            + gamma_final[:, None, None] ** 2
        )
    )
    probabilities = boltz_dist(energies_i, temperature)
    spectrum = np.zeros(
        (len(incident_energies), len(energy_losses), len(pol_type)),
        dtype=float,
    )
    poles = [
        [None for _ in pol_type]
        for _ in incident_energies
    ] if return_poles else None
    retained_in_final_basis = vectors_f.conj().T @ vectors_i
    ntrans = transitions.shape[0]

    polarizations = [
        _rixs_polarization_vectors(
            ntrans, thin, thout, phi,
            incoming_kind, alpha, outgoing_kind, beta, axis,
        )
        for incoming_kind, alpha, outgoing_kind, beta in pol_type
    ]
    for incident_index, omega in enumerate(incident_energies):
        denominator = 1.0 / (
            omega - (energies_n[:, None] - energies_i[None, :])
            + 1j * gamma_core[incident_index]
        )
        for channel, (incoming_pol, outgoing_pol) in enumerate(polarizations):
            absorption_amplitude = np.einsum(
                'k,kni->ni', incoming_pol, absorption, optimize=True
            )
            emission_amplitude = np.einsum(
                'k,kfn->fn', outgoing_pol.conj(), emission, optimize=True
            )
            amplitudes = emission_amplitude @ (
                absorption_amplitude * denominator
            )
            if skip_gs:
                amplitudes -= retained_in_final_basis @ (
                    retained_in_final_basis.conj().T @ amplitudes
                )

            spectrum[incident_index, :, channel] = np.sum(
                final_lorentzian * np.abs(amplitudes)[None, :, :] ** 2
                * probabilities[None, None, :],
                axis=(1, 2),
            )
            if return_poles:
                final_vectors = vectors_f @ amplitudes
                poles[incident_index][channel] = _pole_dict([
                    _pole_record(
                        hamiltonian_i,
                        final_vectors[:, state],
                        energies_i[state],
                    )
                    for state in range(len(energies_i))
                ])

    return (spectrum, poles) if return_poles else spectrum
