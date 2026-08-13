"""Independent tiny-matrix references for SciPy operator unit tests."""

from itertools import combinations

import numpy as np

from edrixs.fock_basis import FockBasis


def fixed_particle_basis(norbs: int, nocc: int) -> FockBasis:
    """Build a small fixed-particle basis used by independent checks."""
    states = []
    for occupied in combinations(range(norbs), nocc):
        state = sum(1 << (norbs - 1 - orbital) for orbital in occupied)
        states.append(state)
    return FockBasis(states, norbs)


def annihilation_operators(norbs: int) -> list[np.ndarray]:
    """Build full-space Jordan-Wigner annihilation matrices."""
    identity = np.eye(2, dtype=complex)
    parity = np.diag([1.0, -1.0]).astype(complex)
    annihilate = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)

    operators = []
    for orbital in range(norbs):
        factors = (
            [parity] * orbital
            + [annihilate]
            + [identity] * (norbs - orbital - 1)
        )
        operator = factors[0]
        for factor in factors[1:]:
            operator = np.kron(operator, factor)
        operators.append(operator)
    return operators


def restrict_operator(full, left_basis: FockBasis, right_basis: FockBasis):
    """Restrict a full Fock-space matrix to selected particle sectors."""
    return full[np.ix_(left_basis.basis_int, right_basis.basis_int)]


def one_body_oracle(emat, left_basis: FockBasis, right_basis=None):
    """Construct an independent one-body Fock-space operator."""
    if right_basis is None:
        right_basis = left_basis
    annihilators = annihilation_operators(left_basis.norbs)
    full = np.zeros_like(annihilators[0])
    for iorb in range(left_basis.norbs):
        for jorb in range(left_basis.norbs):
            full += (
                emat[iorb, jorb]
                * annihilators[iorb].conj().T
                @ annihilators[jorb]
            )
    return restrict_operator(full, left_basis, right_basis)


def two_body_oracle(umat, left_basis: FockBasis, right_basis=None):
    """Construct an independent two-body Fock-space operator."""
    if right_basis is None:
        right_basis = left_basis
    annihilators = annihilation_operators(left_basis.norbs)
    full = np.zeros_like(annihilators[0])
    norbs = left_basis.norbs
    for lorb in range(norbs):
        for korb in range(norbs):
            for jorb in range(norbs):
                for iorb in range(norbs):
                    if umat[lorb, korb, jorb, iorb] != 0:
                        full += (
                            umat[lorb, korb, jorb, iorb]
                            * annihilators[lorb].conj().T
                            @ annihilators[korb].conj().T
                            @ annihilators[jorb]
                            @ annihilators[iorb]
                        )
    return restrict_operator(full, left_basis, right_basis)
