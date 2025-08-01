#!/usr/bin/env python

__all__ = ['fock_bin', 'get_fock_bin_by_N', 'get_fock_half_N',
           'get_fock_full_N', 'get_fock_basis_by_NLz', 'get_fock_basis_by_NSz',
           'get_fock_basis_by_NJz', 'get_fock_basis_by_N_abelian',
           'get_fock_basis_by_N_LzSz', 'write_fock_dec_by_N']

import numpy as np
import itertools


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

    res = []
    for bits in itertools.combinations(list(range(n)), k):
        s = [0] * n
        for bit in bits:
            s[bit] = 1
        res.append(s)
    return res


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

    n = len(args)

    if n % 2 != 0:
        print("Error: number of arguments is not even")
        return

    if n == 2:
        return fock_bin(args[0], args[1])
    else:
        result = []
        res1 = fock_bin(args[0], args[1])
        res2 = get_fock_bin_by_N(*args[2:])
        for ifock in res2:
            for jfock in res1:
                result.append(jfock + ifock)
        return result


def get_fock_half_N(N):
    res = [[] for i in range(N + 1)]
    for i in range(2**N):
        occu = bin(i).count('1')
        res[occu].append(i)
    return res


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
    >>> edrixs.fock_bin(4,2)
    [[1, 1, 0, 0],
     [1, 0, 1, 0],
     [0, 1, 1, 0],
     [1, 0, 0, 1],
     [0, 1, 0, 1],
     [0, 0, 1, 1]]

    >>> import edrixs
    >>> edrixs.get_fock_full_N(4,2)
    [3, 5, 6, 9, 10, 12]

    """

    res = []
    half_N = get_fock_half_N(norb // 2)
    for m in range(norb // 2 + 1):
        n = N - m
        if n >= 0 and n <= norb // 2:
            res.extend([i * 2**(norb // 2) + j for i in half_N[m] for j in half_N[n]])
    return res


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
    {
     -2: [3],
     -1: [5, 6, 9, 10],
      0: [12, 17, 18, 33, 34],
      1: [20, 36, 24, 40],
      2: [48]
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, lz_list)
    return res


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
    {
     -2: [10, 34, 40],
     -1: [],
      0: [3, 6, 9, 12, 18, 33, 36, 24, 48],
      1: [],
      2: [5, 17, 20]
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, sz_list)
    return res


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
    {
     -6: [],
     -5: [],
     -4: [5, 12],
     -3: [],
     -2: [6, 9, 20],
     -1: [],
      0: [3, 10, 17, 36, 24],
      1: [],
      2: [18, 33, 40],
      3: [],
      4: [34, 48],
      5: [],
      6: []
    }
    """

    res = get_fock_basis_by_N_abelian(norb, N, jz_list)
    return res


def get_fock_basis_by_N_abelian(norb, N, a_list):
    """
    Get decimal digitals to represent Fock states, use some Abelian good quantum number.

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
        A dictionary containing the decimal digitals, the key is good quantum numbers,
        the value is a list of int.
    """

    result = get_fock_full_N(norb, N)
    min_a, max_a = min(a_list) * N, max(a_list) * N
    basis = {}
    for i in range(min_a, max_a + 1):
        basis[i] = []
    for n in result:
        a = sum([a_list[i] for i in range(0, n.bit_length()) if (n >> i & 1)])
        basis[a].append(n)
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
        A dictionary containing the decimal digitals, the key is a tuple containing good quantum
        numbers ( :math:`l_{z}`, :math:`s_{z}`), the value is a list of int.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.get_fock_basis_by_N_LzSz(6, 2, [-1, -1, 0, 0, 1, 1], [1, -1, 1, -1, 1, -1])
    {
     (-2, -2): [],
     (-2, -1): [],
      (-2, 0): [3],
      (-2, 1): [],
      (-2, 2): [],
     (-1, -2): [10],
     (-1, -1): [],
      (-1, 0): [6, 9],
      (-1, 1): [],
      (-1, 2): [5],
      (0, -2): [34],
      (0, -1): [],
       (0, 0): [12, 18, 33],
       (0, 1): [],
       (0, 2): [17],
      (1, -2): [40],
      (1, -1): [],
       (1, 0): [36, 24],
       (1, 1): [],
       (1, 2): [20],
      (2, -2): [],
      (2, -1): [],
       (2, 0): [48],
       (2, 1): [],
       (2, 2): []
    }
    """
    result = get_fock_full_N(norb, N)
    min_Lz, max_Lz = min(lz_list) * N, max(lz_list) * N
    min_Sz, max_Sz = min(sz_list) * N, max(sz_list) * N
    basis = {}
    for i in range(min_Lz, max_Lz + 1):
        for j in range(min_Sz, max_Sz + 1):
            basis[(i, j)] = []
    for n in result:
        Lz, Sz = np.sum([[lz_list[i], sz_list[i]] for i in range(0, n.bit_length())
                         if (n >> i & 1)], axis=0)
        basis[(Lz, Sz)].append(n)
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
        The dimension of the Hilbert space

    Examples
    --------
    >>> import edrixs
    >>> edrixs.write_fock_dec_by_N(4, 2, 'fock_i.in')
    file fock_i.in looks like
    15
    3
    5
    6
    9
    10
    12
    17
    18
    20
    24
    33
    34
    36
    40
    48

    where, the first line is the total numer of Fock states,
    and the following lines are the Fock states in decimal form.
    """

    res = get_fock_full_N(N, r)
    res.sort()
    ndim = len(res)
    f = open(fname, 'w')
    print(ndim, file=f)
    for item in res:
        print(item, file=f)
    f.close()
    return ndim
