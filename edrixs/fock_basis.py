#!/usr/bin/env python
"""Fock-basis construction and integer-encoded basis objects."""

from __future__ import annotations

from dataclasses import dataclass, field
import itertools
from math import comb, prod

import numpy as np

__all__ = [
    'FockBasisSpec', 'FockBasis', 'FockBinByN',
    'build_fock_basis', 'get_fock_basis_int', 'get_fock_basis_combinadic',
    'fock_bin', 'get_fock_bin_by_N', 'get_fock_half_N',
    'get_fock_full_N', 'get_fock_basis_by_NLz', 'get_fock_basis_by_NSz',
    'get_fock_basis_by_NJz', 'get_fock_basis_by_N_abelian',
    'get_fock_basis_by_N_LzSz', 'write_fock_dec_by_N',
]


@dataclass(frozen=True, slots=True)
class FockBasisSpec:
    """Compact description of fixed occupancies in one or more orbital shells."""

    shapes: tuple[tuple[int, int], ...]

    def __post_init__(self):
        shapes = tuple((int(norb), int(nocc)) for norb, nocc in self.shapes)
        object.__setattr__(self, 'shapes', shapes)

    @classmethod
    def from_args(cls, *args):
        """Build a specification from ``(norbs, nocc)`` pairs."""
        if len(args) % 2 != 0:
            raise ValueError("number of basis arguments must be even")
        return cls(tuple(zip(args[0::2], args[1::2])))

    @property
    def norbs(self):
        """Total number of spin-orbitals."""
        return sum(norb for norb, _ in self.shapes)

    def __len__(self):
        return prod(comb(norb, nocc) for norb, nocc in self.shapes)


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
    spec : FockBasisSpec or None, optional
        Structured sector description when the basis was generated from fixed
        shell occupancies. Arbitrary explicit bases may leave this as ``None``.
    """

    def __init__(self, basis_int, norbs, spec=None):
        self.basis_int = [int(value) for value in basis_int]
        self.norbs = int(norbs)
        self.spec = spec
        self.lookup = {value: index for index, value in enumerate(self.basis_int)}

    def __len__(self):
        return len(self.basis_int)

    def encode(self, value):
        """Return the basis position of an integer-encoded Fock state."""
        return self.lookup[int(value)]

    def decode(self, position):
        """Return the integer-encoded Fock state at ``position``."""
        return self.basis_int[position]


@dataclass(slots=True)
class FockBinByN:
    """Implicit combinadic basis for a complete fixed-occupancy sector."""

    shapes: tuple[tuple[int, int], ...]
    spec: FockBasisSpec = field(init=False)
    shell_norbs: tuple[int, ...] = field(init=False)
    norbs: int = field(init=False)
    noccus: tuple[int, ...] = field(init=False)
    sizes: tuple[int, ...] = field(init=False)
    num_subspaces: int = field(init=False)
    num_orbitals: int = field(init=False)
    dim: int = field(init=False)
    offsets: tuple[int, ...] = field(init=False)
    strides: tuple[int, ...] = field(init=False)
    min_decode: int = field(init=False)
    max_decode: int = field(init=False)

    def __post_init__(self):
        self.spec = FockBasisSpec(tuple(self.shapes))
        self.shapes = self.spec.shapes
        self.shell_norbs = tuple(norb for norb, _ in self.shapes)
        self.norbs = self.spec.norbs
        self.noccus = tuple(nocc for _, nocc in self.shapes)
        self.sizes = tuple(comb(norb, nocc) for norb, nocc in self.shapes)
        self.num_subspaces = len(self.shapes)
        self.num_orbitals = self.norbs
        self.dim = prod(self.sizes)

        running = self.norbs
        offsets = []
        for norb in self.shell_norbs:
            running -= norb
            offsets.append(running)
        self.offsets = tuple(offsets)

        # The first shell varies fastest, matching get_fock_bin_by_N.
        stride = 1
        strides = []
        for size in self.sizes:
            strides.append(stride)
            stride *= size
        self.strides = tuple(strides)
        self.min_decode, self.max_decode = _min_max_decode(self.shapes)

    @classmethod
    def from_spec(cls, spec):
        """Build an implicit basis from a :class:`FockBasisSpec`."""
        return cls(spec.shapes)

    def __len__(self):
        return self.dim

    def encode(self, state):
        """Return the canonical EDRIXS index of ``state``."""
        index = _encode_combinadic(
            int(state), self.shell_norbs, self.noccus, self.offsets,
            self.sizes, self.strides,
        )
        if index < 0:
            raise KeyError(int(state))
        return index

    def decode(self, index):
        """Return the integer-encoded state at the canonical EDRIXS index."""
        if index < 0:
            index += self.dim
        if index < 0 or index >= self.dim:
            raise IndexError(index)
        return _decode_combinadic(
            int(index), self.shell_norbs, self.noccus, self.offsets, self.sizes
        )

    def jit_args(self):
        """Return compact numeric metadata for an explicitly requested JIT path."""
        return (
            np.asarray(self.shell_norbs, dtype=np.int64),
            np.asarray(self.noccus, dtype=np.int64),
            np.asarray(self.offsets, dtype=np.int64),
            np.asarray(self.sizes, dtype=np.int64),
            np.asarray(self.strides, dtype=np.int64),
            np.uint64(self.min_decode),
            np.uint64(self.max_decode),
        )


def _hash_decoder(rank, norb, nocc):
    """Decode a conventional colex combinadic rank inside one shell."""
    if nocc == 0:
        return 0
    state = 0
    j = nocc
    i = norb - 1
    c = comb(i, j)
    while i >= 0 and j > 0:
        if c <= rank:
            state |= 1 << i
            rank -= c
            old_i, old_j = i, j
            i -= 1
            j -= 1
            if j == 0 or i < 0:
                break
            c = (c * old_j) // old_i
        else:
            if i == 0:
                break
            c = (c * (i - j)) // i
            i -= 1
    return state


def _hash_encoder(state, norb):
    """Encode one shell into its conventional colex combinadic rank."""
    state = int(state) & ((1 << norb) - 1)
    rank = 0
    k = 1
    while state:
        lsb = state & -state
        pos = lsb.bit_length() - 1
        rank += comb(pos, k)
        k += 1
        state ^= lsb
    return rank


def _decode_combinadic(index, norbs, noccus, offsets, sizes):
    """Decode using the historical EDRIXS shell and state ordering."""
    state = 0
    for norb, nocc, offset, size in zip(norbs, noccus, offsets, sizes):
        shell_rank = index % size
        index //= size
        colex_rank = size - 1 - shell_rank
        state |= _hash_decoder(colex_rank, norb, nocc) << offset
    return state


def _encode_combinadic(state, norbs, noccus, offsets, sizes, strides):
    """Encode using the historical EDRIXS shell and state ordering."""
    index = 0
    for norb, nocc, offset, size, stride in zip(
            norbs, noccus, offsets, sizes, strides):
        mask = (1 << norb) - 1
        shell_state = (state >> offset) & mask
        if shell_state.bit_count() != nocc:
            return -1
        shell_rank = size - 1 - _hash_encoder(shell_state, norb)
        index += shell_rank * stride
    return index


def _min_max_decode(shapes):
    total_norb = sum(norb for norb, _ in shapes)
    state_min = 0
    state_max = 0
    running = total_norb
    for norb, nocc in shapes:
        running -= norb
        shell_min = (1 << nocc) - 1
        shell_max = shell_min << (norb - nocc)
        state_min |= shell_min << running
        state_max |= shell_max << running
    return state_min, state_max


def get_fock_basis_int(*args):
    """
    Build an explicit :class:`FockBasis` for fixed shell occupancies.

    The state ordering is the historical EDRIXS ordering produced by
    :func:`get_fock_bin_by_N`.
    """
    if len(args) % 2 != 0:
        print("Error: number of arguments is not even")
        return None
    spec = FockBasisSpec.from_args(*args)
    basis_binary = get_fock_bin_by_N(*args)
    basis_int = np.asarray(
        [int(''.join(map(str, row)), 2) for row in basis_binary],
        dtype=object,
    )
    return FockBasis(basis_int, spec.norbs, spec=spec)


def get_fock_basis_combinadic(*args):
    """Build an implicit combinadic basis for fixed shell occupancies."""
    if len(args) % 2 != 0:
        print("Error: number of arguments is not even")
        return None
    return FockBinByN.from_spec(FockBasisSpec.from_args(*args))


def build_fock_basis(basis, method='combinadic'):
    """Realize a compact basis specification with the requested representation."""
    if isinstance(basis, (FockBasis, FockBinByN)):
        return basis
    if not isinstance(basis, FockBasisSpec):
        raise TypeError("basis must be a FockBasisSpec or a realized Fock basis")

    if method == 'combinadic':
        return FockBinByN.from_spec(basis)
    if method == 'explicit':
        args = tuple(value for shape in basis.shapes for value in shape)
        return get_fock_basis_int(*args)
    raise ValueError("basis method must be 'combinadic' or 'explicit'")


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
