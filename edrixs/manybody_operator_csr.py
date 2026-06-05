import numpy as np
import scipy.sparse as sp

try:
    from .fock_basis import get_fock_bin_by_N
except ImportError:  # pragma: no cover
    from fock_basis import get_fock_bin_by_N


class FockBasis:
    """
    Integer-encoded Fock basis.

    The occupation pattern is stored as an integer bit string.  Orbital 0 is
    the most significant bit, matching the ordering used by get_fock_bin_by_N.

    Parameters
    ----------
    basis_int : sequence of int
        Integer-encoded occupation patterns.

    norbs : int
        Number of spin-orbitals.

    Examples
    --------
    If the occupation vector is [1, 1, 0, 0], the encoded integer is
    int('1100', 2) = 12.
    """

    def __init__(self, basis_int, norbs):
        self.basis_int = [int(v) for v in basis_int]
        self.norbs = int(norbs)
        self.lookup = {value: i for i, value in enumerate(self.basis_int)}

    def __len__(self):
        return len(self.basis_int)

    def encode(self, value):
        """Return the basis position of an integer-encoded Fock state."""
        return self.lookup[int(value)]

    def decode(self, position):
        """Return the integer-encoded Fock state at a basis position."""
        return self.basis_int[position]


def get_fock_basis_int(*args):
    """
    Get an integer-encoded Fock basis for a multi-shell system.

    Parameters
    ----------
    args : ints
        Pairs of ``number_of_orbitals, occupancy`` for each shell.

    Returns
    -------
    basis : FockBasis
        Integer-encoded Fock basis with ``encode`` and ``decode`` methods.
    """
    basis_10 = get_fock_bin_by_N(*args)
    basis_int = np.array(
        [int(''.join(map(str, row)), 2) for row in basis_10],
        dtype=object
    )
    norbs = sum(args[0::2])
    return FockBasis(basis_int, norbs)


def _count_occupied_before(state, orb, norbs):
    """Count occupied orbitals with index smaller than orb."""
    return (int(state) >> (norbs - orb)).bit_count()


def two_fermion_csr(emat, lb, rb=None, tol=1E-10):
    """
    Build matrix form of a two-fermion operator in the given Fock basis.

    The operator convention is

    .. math::

        <F_l|\\sum_{ij} E_{ij} f_i^\\dagger f_j |F_r>.

    Parameters
    ----------
    emat : 2d complex array
        One-body orbital matrix.
    lb : FockBasis
        Left Fock basis :math:`<F_l|`.
    rb : FockBasis, optional
        Right Fock basis :math:`|F_r>`.  If None, rb = lb.
    tol : float, optional
        Ignore elements of emat with absolute value <= tol.

    Returns
    -------
    hmat : scipy.sparse.csr_matrix
        Sparse many-body representation of the one-body operator.
    """
    if rb is None:
        rb = lb

    nr = len(rb)
    nl = len(lb)
    norbs = lb.norbs

    if rb.norbs != norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    a1, a2 = np.nonzero(np.abs(emat) > tol)
    nonzero = np.stack((a1, a2), axis=-1)

    rows = []
    cols = []
    data = []

    for iorb, jorb in nonzero:
        bit_j = 1 << (norbs - 1 - jorb)
        bit_i = 1 << (norbs - 1 - iorb)

        for icfg in range(nr):
            state = rb.decode(icfg)

            if not state & bit_j:
                continue
            s1 = (-1) ** _count_occupied_before(state, jorb, norbs)
            state ^= bit_j

            if state & bit_i:
                continue
            s2 = (-1) ** _count_occupied_before(state, iorb, norbs)
            state |= bit_i

            try:
                jcfg = lb.encode(state)
            except KeyError:
                continue

            rows.append(jcfg)
            cols.append(icfg)
            data.append(emat[iorb, jorb] * s1 * s2)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(nl, nr),
        dtype=np.complex128
    ).tocsr()


def four_fermion_csr(umat, lb, rb=None, tol=1E-10):
    """
    Build matrix form of a four-fermion operator in the given Fock basis.

    The tensor convention is

    .. math::

        <F_l|\\sum_{lkji} U_{lkji}
        f_l^\\dagger f_k^\\dagger f_j f_i |F_r>.

    Parameters
    ----------
    umat : 4d complex array
        Coulomb tensor with entries ``umat[lorb, korb, jorb, iorb]``.
    lb : FockBasis
        Left Fock basis :math:`<F_l|`.
    rb : FockBasis, optional
        Right Fock basis :math:`|F_r>`.  If None, rb = lb.
    tol : float, optional
        Ignore elements of umat with absolute value <= tol.

    Returns
    -------
    hmat : scipy.sparse.csr_matrix
        Sparse many-body representation of the two-body operator.
    """
    if rb is None:
        rb = lb

    nr = len(rb)
    nl = len(lb)
    norbs = lb.norbs

    if rb.norbs != norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    a1, a2, a3, a4 = np.nonzero(np.abs(umat) > tol)
    nonzero = np.stack((a1, a2, a3, a4), axis=-1)

    rows = []
    cols = []
    data = []

    for lorb, korb, jorb, iorb in nonzero:
        if iorb == jorb or korb == lorb:
            continue

        bit_i = 1 << (norbs - 1 - iorb)
        bit_j = 1 << (norbs - 1 - jorb)
        bit_k = 1 << (norbs - 1 - korb)
        bit_l = 1 << (norbs - 1 - lorb)

        for icfg in range(nr):
            state = rb.decode(icfg)

            if not state & bit_i:
                continue
            s1 = (-1) ** _count_occupied_before(state, iorb, norbs)
            state ^= bit_i

            if not state & bit_j:
                continue
            s2 = (-1) ** _count_occupied_before(state, jorb, norbs)
            state ^= bit_j

            if state & bit_k:
                continue
            s3 = (-1) ** _count_occupied_before(state, korb, norbs)
            state |= bit_k

            if state & bit_l:
                continue
            s4 = (-1) ** _count_occupied_before(state, lorb, norbs)
            state |= bit_l

            try:
                jcfg = lb.encode(state)
            except KeyError:
                continue

            rows.append(jcfg)
            cols.append(icfg)
            data.append(umat[lorb, korb, jorb, iorb] * s1 * s2 * s3 * s4)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(nl, nr),
        dtype=np.complex128
    ).tocsr()


def four_fermion_csr_auto(umat, basis, rb=None, tol=1E-10):
    """Dispatch four-fermion construction based on dense or sparse U format."""
    if sp.issparse(umat):
        return _four_fermion_csr_from_sparse_umat(umat, basis, rb=rb, tol=tol)

    return four_fermion_csr(umat, basis, rb=rb, tol=tol)


def _four_fermion_csr_from_sparse_umat(umat, lb, rb=None, tol=1E-10):
    """
    Build a four-fermion many-body operator from sparse flattened U.

    The sparse U matrix must have shape (norbs*norbs, norbs*norbs), with

        row = lorb * norbs + korb
        col = jorb * norbs + iorb

    representing tensor entries umat[lorb, korb, jorb, iorb].
    """
    if rb is None:
        rb = lb

    nr = len(rb)
    nl = len(lb)
    norbs = lb.norbs

    if rb.norbs != norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    expected_shape = (norbs * norbs, norbs * norbs)
    if umat.shape != expected_shape:
        raise ValueError(
            "sparse umat has shape {}, expected {}".format(
                umat.shape, expected_shape
            )
        )

    umat = umat.tocoo()

    rows = []
    cols = []
    data = []

    for row, col, val in zip(umat.row, umat.col, umat.data):
        if abs(val) <= tol:
            continue

        lorb = row // norbs
        korb = row % norbs
        jorb = col // norbs
        iorb = col % norbs

        if iorb == jorb or korb == lorb:
            continue

        bit_i = 1 << (norbs - 1 - iorb)
        bit_j = 1 << (norbs - 1 - jorb)
        bit_k = 1 << (norbs - 1 - korb)
        bit_l = 1 << (norbs - 1 - lorb)

        for icfg in range(nr):
            state = rb.decode(icfg)

            if not state & bit_i:
                continue
            s1 = (-1) ** (int(state) >> (norbs - iorb)).bit_count()
            state ^= bit_i

            if not state & bit_j:
                continue
            s2 = (-1) ** (int(state) >> (norbs - jorb)).bit_count()
            state ^= bit_j

            if state & bit_k:
                continue
            s3 = (-1) ** (int(state) >> (norbs - korb)).bit_count()
            state |= bit_k

            if state & bit_l:
                continue
            s4 = (-1) ** (int(state) >> (norbs - lorb)).bit_count()
            state |= bit_l

            try:
                jcfg = lb.encode(state)
            except KeyError:
                continue

            rows.append(jcfg)
            cols.append(icfg)
            data.append(val * s1 * s2 * s3 * s4)

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(nl, nr),
        dtype=np.complex128
    ).tocsr()
