"""Small independent dense oracles used only by the tests."""

from itertools import combinations

import numpy as np

from edrixs.manybody_operator_csr import FockBasis


def fixed_particle_basis(norbs: int, nocc: int) -> FockBasis:
    states = []
    for occupied in combinations(range(norbs), nocc):
        state = sum(1 << (norbs - 1 - orb) for orb in occupied)
        states.append(state)
    return FockBasis(states, norbs)


def annihilation_operators(norbs: int) -> list[np.ndarray]:
    """Jordan-Wigner annihilation matrices; orbital 0 is the MSB."""
    ident = np.eye(2, dtype=complex)
    parity = np.diag([1.0, -1.0]).astype(complex)
    annihilate = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)

    operators = []
    for orb in range(norbs):
        factors = [parity] * orb + [annihilate] + [ident] * (norbs - orb - 1)
        op = factors[0]
        for factor in factors[1:]:
            op = np.kron(op, factor)
        operators.append(op)
    return operators


def restrict_operator(full, left_basis: FockBasis, right_basis: FockBasis):
    return full[np.ix_(left_basis.basis_int, right_basis.basis_int)]


def one_body_oracle(emat, left_basis: FockBasis, right_basis=None):
    if right_basis is None:
        right_basis = left_basis
    c = annihilation_operators(left_basis.norbs)
    full = np.zeros_like(c[0])
    for i in range(left_basis.norbs):
        for j in range(left_basis.norbs):
            full += emat[i, j] * c[i].conj().T @ c[j]
    return restrict_operator(full, left_basis, right_basis)


def two_body_oracle(umat, left_basis: FockBasis, right_basis=None):
    if right_basis is None:
        right_basis = left_basis
    c = annihilation_operators(left_basis.norbs)
    full = np.zeros_like(c[0])
    n = left_basis.norbs
    for lorb in range(n):
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    if umat[lorb, k, j, i] != 0:
                        full += (
                            umat[lorb, k, j, i]
                            * c[lorb].conj().T
                            @ c[k].conj().T
                            @ c[j]
                            @ c[i]
                        )
    return restrict_operator(full, left_basis, right_basis)
