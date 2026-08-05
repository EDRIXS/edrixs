#!/usr/bin/env python
"""Fock-basis construction and integer-encoded basis objects."""

from __future__ import annotations

import itertools
import numpy as np

__all__ = [
    'FockBasis', 'get_fock_basis_int',
    'fock_bin', 'get_fock_bin_by_N', 'get_fock_half_N',
    'get_fock_full_N', 'get_fock_basis_by_NLz', 'get_fock_basis_by_NSz',
    'get_fock_basis_by_NJz', 'get_fock_basis_by_N_abelian',
    'get_fock_basis_by_N_LzSz', 'write_fock_dec_by_N',
]


class FockBasis:
    """Integer-encoded Fock basis.

    Orbital zero is stored as the most significant bit, matching
    :func:`get_fock_bin_by_N` and the many-body construction backends.
    """

    def __init__(self, basis_int, norbs):
        self.basis_int = [int(value) for value in basis_int]
        self.norbs = int(norbs)
        self.lookup = {value: index for index, value in enumerate(self.basis_int)}

    def __len__(self):
        return len(self.basis_int)

    def encode(self, value):
        """Return the basis position of an integer-encoded Fock state."""
        return self.lookup[int(value)]

    def decode(self, position):
        """Return the integer-encoded Fock state at ``position``."""
        return self.basis_int[position]


def get_fock_basis_int(*args):
    """Return an integer-encoded :class:`FockBasis` for shell pairs.

    Parameters are ``(number_of_orbitals, occupancy)`` pairs, with the same
    ordering convention as :func:`get_fock_bin_by_N`.
    """
    basis_binary = get_fock_bin_by_N(*args)
    if basis_binary is None:
        return None
    basis_int = np.asarray(
        [int(''.join(map(str, row)), 2) for row in basis_binary],
        dtype=object,
    )
    return FockBasis(basis_int, sum(args[0::2]))


def fock_bin(n, k):
    """Return all length-``n`` binary occupations containing ``k`` ones."""
    if n == 0:
        return [[0]]

    result = []
    for occupied in itertools.combinations(range(n), k):
        state = [0] * n
        for orbital in occupied:
            state[orbital] = 1
        result.append(state)
    return result


def get_fock_bin_by_N(*args):
    """Return binary Fock states for one or more fixed-occupancy shells."""
    narg = len(args)
    if narg % 2 != 0:
        print("Error: number of arguments is not even")
        return None

    if narg == 2:
        return fock_bin(args[0], args[1])

    result = []
    first_shell = fock_bin(args[0], args[1])
    remaining_shells = get_fock_bin_by_N(*args[2:])
    for remaining in remaining_shells:
        for first in first_shell:
            result.append(first + remaining)
    return result


def get_fock_half_N(N):
    result = [[] for _ in range(N + 1)]
    for state in range(2**N):
        result[bin(state).count('1')].append(state)
    return result


def get_fock_full_N(norb, N):
    """Return decimal Fock states with total occupancy ``N``."""
    result = []
    half_basis = get_fock_half_N(norb // 2)
    for left_occupancy in range(norb // 2 + 1):
        right_occupancy = N - left_occupancy
        if 0 <= right_occupancy <= norb // 2:
            result.extend(
                left * 2**(norb // 2) + right
                for left in half_basis[left_occupancy]
                for right in half_basis[right_occupancy]
            )
    return result


def get_fock_basis_by_NLz(norb, N, lz_list):
    return get_fock_basis_by_N_abelian(norb, N, lz_list)


def get_fock_basis_by_NSz(norb, N, sz_list):
    return get_fock_basis_by_N_abelian(norb, N, sz_list)


def get_fock_basis_by_NJz(norb, N, jz_list):
    return get_fock_basis_by_N_abelian(norb, N, jz_list)


def get_fock_basis_by_N_abelian(norb, N, a_list):
    """Group fixed-occupancy Fock states by an Abelian quantum number."""
    result = get_fock_full_N(norb, N)
    minimum, maximum = min(a_list) * N, max(a_list) * N
    basis = {value: [] for value in range(minimum, maximum + 1)}
    for state in result:
        quantum_number = sum(
            a_list[index]
            for index in range(state.bit_length())
            if state >> index & 1
        )
        basis[quantum_number].append(state)
    return basis


def get_fock_basis_by_N_LzSz(norb, N, lz_list, sz_list):
    """Group fixed-occupancy Fock states by ``(Lz, Sz)``."""
    result = get_fock_full_N(norb, N)
    min_lz, max_lz = min(lz_list) * N, max(lz_list) * N
    min_sz, max_sz = min(sz_list) * N, max(sz_list) * N
    basis = {
        (lz, sz): []
        for lz in range(min_lz, max_lz + 1)
        for sz in range(min_sz, max_sz + 1)
    }
    for state in result:
        occupied = [
            index for index in range(state.bit_length())
            if state >> index & 1
        ]
        lz = sum(lz_list[index] for index in occupied)
        sz = sum(sz_list[index] for index in occupied)
        basis[(lz, sz)].append(state)
    return basis


def write_fock_dec_by_N(N, r, fname='fock_i.in'):
    """Write sorted decimal fixed-occupancy Fock states to ``fname``."""
    result = get_fock_full_N(N, r)
    result.sort()
    with open(fname, 'w') as stream:
        print(len(result), file=stream)
        for item in result:
            print(item, file=stream)
    return len(result)
