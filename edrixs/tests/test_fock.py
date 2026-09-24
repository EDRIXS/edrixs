import edrixs


def test_get_fock_full_N():
    """Reproduce known behavior of get_fock_full_N"""
    assert edrixs.get_fock_full_N(4, 2) == [3, 5, 6, 9, 10, 12]


def get_fock_half_N(N):
    res = [[] for i in range(N + 1)]
    for i in range(2**N):
        occu = bin(i).count('1')
        res[occu].append(i)
    return res


def get_fock_full_N_old(norb, N):
    """
    The slower function only works for even norb.
    Use this to test the revised function.

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
    """

    res = []
    half_N = get_fock_half_N(norb // 2)
    for m in range(norb // 2 + 1):
        n = N - m
        if n >= 0 and n <= norb // 2:
            res.extend([i * 2**(norb // 2) + j
                        for i in half_N[m] for j in half_N[n]])
    return res


def test_old_new():
    for norb in [6, 8, 10]:  # only test even norb
        N = norb - 3
        assert (sorted(get_fock_full_N_old(norb, N)) ==
                sorted(edrixs.get_fock_full_N(norb, N)))


def binary_from_full(norb, N):
    """Helper function, which generates binary
    from the fock_bin function."""
    return [int(''.join(str(b) for b in binary), 2)
            for binary in edrixs.fock_bin(norb, N)]


def test_fast_fock_compare():
    """Confirm that fast version of the function behaves
    the same as the more easily interpretable slow version."""
    for norb in [10, 15, 20]:
        assert (sorted(edrixs.get_fock_full_N(norb, norb - 5))
                == sorted(binary_from_full(norb, norb - 5)))
