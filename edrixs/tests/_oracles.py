"""Independent tiny-matrix references for sparse operator unit tests.

These helpers build one- and two-body many-electron operators without using
``ops(..., backend="scipy")``.  They let the unit tests check the stage that
turns orbital-space terms from ``setup_*`` into sparse Fock-space operators
before those operators reach the ED, XAS, or RIXS solvers.
"""

from itertools import combinations

import numpy as np

from edrixs.manybody_operator_csr import FockBasis


def fixed_particle_basis(norbs: int, nocc: int) -> FockBasis:
    """Build a small fixed-particle basis used by independent operator checks."""
    states = []
    for occupied in combinations(range(norbs), nocc):
        state = sum(1 << (norbs - 1 - orb) for orb in occupied)
        states.append(state)
    return FockBasis(states, norbs)


def annihilation_operators(norbs: int) -> list[np.ndarray]:
    """Build full-space Jordan-Wigner annihilation matrices for each orbital."""
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
    """Restrict a full Fock-space matrix to the particle sectors under test."""
    return full[np.ix_(left_basis.basis_int, right_basis.basis_int)]


def one_body_oracle(emat, left_basis: FockBasis, right_basis=None):
    """Construct an independent one-body Fock-space operator for comparison."""
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
    """Construct an independent Coulomb operator for sparse-path comparison."""
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
