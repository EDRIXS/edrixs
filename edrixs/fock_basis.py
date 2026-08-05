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
    """
    Integer-encoded Fock basis with constant-time state lookup.

    Orbital zero is stored as the most significant bit, matching
    :func:`get_fock_bin_by_N` and the many-body operator backends.

    Parameters
    ----------
    basis_int : sequence of int
        Integer encodings of the basis states in matrix-index order.
    norbs : int
        Number of spin-orbitals represented by each state.
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
    """
    Build an integer-encoded :class:`FockBasis` for fixed shell occupancies.

    Parameters
    ----------
    args : ints
        ``(number_of_orbitals, occupancy)`` pairs, with the same shell ordering
        convention as :func:`get_fock_bin_by_N`.

    Returns
    -------
    FockBasis or None
        Integer-encoded basis, or ``None`` when an odd number of arguments is
        supplied.
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
    """
    Return all the possible :math:`n`-length binary
    where :math:`k` of :math:`n` digitals are set to 1.

    Parameters
    ----------
    n: int
        Binary length :math:`n`.
    k: int
        How many digitals are set to be 1.

    Returns
    -------
    res: list of int-lists
        A list of list containing the binary digitals.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.fock_bin(4, 2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 1, 0],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]
    """
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
    """
    Get binary form to represent a Fock state.

    Parameters
    ----------
    args: ints
        args[0]: number of orbitals for 1st-shell,

        args[1]: number of occupancy for 1st-shell,

        args[2]: number of orbitals for 2nd-shell,

        args[3]: number of occupancy for 2nd-shell,

        ...

        args[ :math:`2N-2`]: number of orbitals for :math:`N` th-shell,

        args[ :math:`2N-1`]: number of occupancy for :math:`N` th-shell.

    Returns
    -------
    result: list of int list
        The binary form of Fock states.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_bin_by_N(4, 2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 1, 0],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    >>> edrixs.get_fock_bin_by_N(4, 2, 2, 1)
    [[1, 1, 0, 0, 1, 0],
     [1, 0, 1, 0, 1, 0],
     [1, 0, 0, 1, 1, 0],
     [0, 1, 1, 0, 1, 0],
     [0, 1, 0, 1, 1, 0],
     [0, 0, 1, 1, 1, 0],
     [1, 1, 0, 0, 0, 1],
     [1, 0, 1, 0, 0, 1],
     [1, 0, 0, 1, 0, 1],
     [0, 1, 1, 0, 0, 1],
     [0, 1, 0, 1, 0, 1],
     [0, 0, 1, 1, 0, 1]]
    """
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
    """Group half-space integer states by particle number."""
    result = [[] for _ in range(N + 1)]
    for state in range(2**N):
        result[bin(state).count('1')].append(state)
    return result


def get_fock_full_N(norb, N):
    """
    Get the decimal digitals to represent Fock states.

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of occupancy.

    Returns
    -------
    res: list of int
        The decimal digitals to represent Fock states.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.fock_bin(4, 2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [0, 1, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    >>> edrixs.get_fock_full_N(4, 2)
    [3, 5, 6, 9, 10, 12]
    """
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
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - orbital angular momentum :math:`L_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    lz_list: list of int
        Quantum number :math:`l_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good
        quantum numbers :math:`L_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NLz(6, 2, [-1, -1, 0, 0, 1, 1])
    {-2: [3], -1: [5, 6, 9, 10], 0: [12, 17, 18, 33, 34],
     1: [20, 36, 24, 40], 2: [48]}
    """
    return get_fock_basis_by_N_abelian(norb, N, lz_list)


def get_fock_basis_by_NSz(norb, N, sz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - spin angular momentum :math:`S_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    sz_list: list of int
        Quantum number :math:`s_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good quantum
        numbers :math:`S_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NSz(6, 2, [1, -1, 1, -1, 1, -1])
    {-2: [10, 34, 40], -1: [],
     0: [3, 6, 9, 12, 18, 33, 36, 24, 48], 1: [], 2: [5, 17, 20]}
    """
    return get_fock_basis_by_N_abelian(norb, N, sz_list)


def get_fock_basis_by_NJz(norb, N, jz_list):
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - total angular momentum :math:`J_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    jz_list: list of int
        Quantum number :math:`j_{z}` for each orbital.

    Returns
    -------
    res: dict
        A dictionary containing the decimal digitals, the key is good quantum
        numbers :math:`j_{z}`, the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_NJz(6, 2, [-1, 1, -3, -1, 1, 3])
    {-6: [], -5: [], -4: [5, 12], -3: [], -2: [6, 9, 20], -1: [],
     0: [3, 10, 17, 36, 24], 1: [], 2: [18, 33, 40], 3: [],
     4: [34, 48], 5: [], 6: []}
    """
    return get_fock_basis_by_N_abelian(norb, N, jz_list)


def get_fock_basis_by_N_abelian(norb, N, a_list):
    """
    Get decimal digitals to represent Fock states, use some Abelian good
    quantum number.

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    a_list: list of int
        Quantum number of the Abelian symmetry for each orbital.

    Returns
    -------
    basis: dict
        A dictionary containing the decimal digitals, the key is good quantum
        numbers, the value is a list of int.
    """
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
    """
    Get decimal digitals to represent Fock states, use good quantum number:

    - orbital angular momentum :math:`L_{z}`
    - spin angular momentum :math:`S_{z}`

    Parameters
    ----------
    norb: int
        Number of orbitals.
    N: int
        Number of total occupancy.
    lz_list: list of int
        Quantum number :math:`l_{z}` for each orbital.
    sz_list: list of int
        Quantum number :math:`s_{z}` for each orbital.

    Returns
    -------
    basis: dict
        A dictionary containing the decimal digitals, the key is a tuple of
        good quantum numbers (:math:`l_{z}`, :math:`s_{z}`), and the value is
        a list of int.

    Examples
    --------
    >>> import edrixs
    >>> lz = [-1, -1, 0, 0, 1, 1]
    >>> sz = [1, -1, 1, -1, 1, -1]
    >>> edrixs.get_fock_basis_by_N_LzSz(6, 2, lz, sz)[(0, 0)]
    [12, 18, 33]
    """
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
    """
    Get decimal digitals to represent Fock states, sort them by
    ascending order and then write them to file.

    Parameters
    ----------
    N: int
       Number of orbitals.
    r: int
        Number of occuancy.
    fname: string
        File name.

    Returns
    -------
    ndim: int
        The dimension of the Hilbert space.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    6

    The first line of ``fock_i.in`` is the total number of Fock states,
    followed by the states in decimal form: ``3, 5, 6, 9, 10, 12``.
    """
    result = get_fock_full_N(N, r)
    result.sort()
    with open(fname, 'w') as stream:
        print(len(result), file=stream)
        for item in result:
            print(item, file=stream)
    return len(result)
