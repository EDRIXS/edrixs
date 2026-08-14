import numpy as np
from math import comb, prod
from dataclasses import dataclass, field
from typing import List, Tuple

from numba import njit
from petsc4py import PETSc


@dataclass(slots=True)
class FockBinByN:
    shapes: List[Tuple[int, int]]
    norbs: List[int] = field(init=False)
    noccus: List[int] = field(init=False)
    sizes: List[int] = field(init=False)
    num_subspaces: int = field(init=False)
    num_orbitals: int = field(init=False)
    dim: int = field(init=False)
    offsets: Tuple[int, ...] = field(init=False)
    strides: Tuple[int, ...] = field(init=False)
    min_decode: int = field(init=False)
    max_decode: int = field(init=False)

    def __post_init__(self) -> None:
        self.norbs = [int(norb) for (norb, _) in self.shapes]
        self.noccus = [int(nocc) for (_, nocc) in self.shapes]
        self.sizes = [comb(norb, nocc) for (norb, nocc) in self.shapes]
        self.num_subspaces = len(self.norbs)
        self.num_orbitals = sum(self.norbs)
        self.dim = prod(self.sizes)

        running = self.num_orbitals
        offsets = []
        for N in self.norbs:
            running -= N
            offsets.append(running)
        self.offsets = tuple(offsets)

        self.min_decode, self.max_decode = min_max_decode(self.shapes)

        stride = 1
        strides = [0] * self.num_subspaces
        for n in reversed(range(self.num_subspaces)):
            strides[n] = stride
            stride *= self.sizes[n]
        self.strides = tuple(strides)

    def encode(self, b: int) -> int:
        return encode_basis_py(b, self.norbs, self.noccus, self.offsets, self.strides,
                               self.min_decode, self.max_decode)

    def decode(self, index: int) -> int:
        return decode_basis_py(index, self.norbs, self.noccus, self.offsets, self.sizes)

    def jit_args(self):
        """
        Small metadata only. No full state table.
        """
        return (
            np.asarray(self.norbs, dtype=np.int64),
            np.asarray(self.noccus, dtype=np.int64),
            np.asarray(self.offsets, dtype=np.int64),
            np.asarray(self.sizes, dtype=np.int64),
            np.asarray(self.strides, dtype=np.int64),
            np.int64(self.min_decode),
            np.int64(self.max_decode),
        )


def hash_decoder(r: int, N: int, M: int) -> int:
    b = 0
    j = M
    i = N - 1
    c = comb(i, j)
    while i >= 0 and j > 0:
        if c <= r:
            b |= (1 << i)
            r -= c
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
    return b


def hash_encoder(b: int, N: int, M: int) -> int:
    b = int(b) & ((1 << N) - 1)
    r = 0
    k = 1
    while b:
        lsb = b & -b
        pos = 0
        t = lsb
        while t > 1:
            t >>= 1
            pos += 1
        r += comb(pos, k)
        k += 1
        b ^= lsb
    return r


def min_max_decode(shapes):
    totalN = sum(N for N, _ in shapes)
    b_min = 0
    b_max = 0
    running = totalN
    for (N, M) in shapes:
        running -= N
        sub_min = (1 << M) - 1
        sub_max = ((1 << M) - 1) << (N - M)
        b_min |= sub_min << running
        b_max |= sub_max << running
    return b_min, b_max


@njit(cache=True, inline="always")
def comb_jit(n, k):
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n - k:
        k = n - k
    c = 1
    for i in range(k):
        c = (c * (n - i)) // (i + 1)
    return c


@njit(cache=True, inline="always")
def popcount_u64(x):
    c = 0
    while x != np.uint64(0):
        x &= x - np.uint64(1)
        c += 1
    return c


@njit(cache=True, inline="always")
def hash_decoder_jit(r, N, M):
    b = np.uint64(0)
    j = M
    i = N - 1
    c = comb_jit(i, j)

    while i >= 0 and j > 0:
        if c <= r:
            b |= np.uint64(1) << np.uint64(i)
            r -= c
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
    return b


@njit(cache=True, inline="always")
def hash_encoder_jit(b, N, M):
    b = np.uint64(b) & ((np.uint64(1) << np.uint64(N)) - np.uint64(1))
    r = 0
    k = 1
    while b != np.uint64(0):
        lsb = b & (np.uint64(0) - b)
        pos = 0
        t = lsb
        while t > np.uint64(1):
            t >>= np.uint64(1)
            pos += 1
        r += comb_jit(pos, k)
        k += 1
        b ^= lsb
    return r


@njit(cache=True, inline="always")
def bit_at_orb_u64(b, orb, norb):
    bitpos = np.uint64(norb - 1 - orb)
    return (b >> bitpos) & np.uint64(1)


@njit(cache=True, inline="always")
def set_zero_at_orb_u64(b, orb, norb):
    bitpos = np.uint64(norb - 1 - orb)
    return b & ~(np.uint64(1) << bitpos)


@njit(cache=True, inline="always")
def set_one_at_orb_u64(b, orb, norb):
    bitpos = np.uint64(norb - 1 - orb)
    return b | (np.uint64(1) << bitpos)


@njit(cache=True, inline="always")
def sign_count_u64(b, orb, norb):
    bitpos = np.uint64(norb - 1 - orb)
    prefix = b >> (bitpos + np.uint64(1))
    return 1 if (popcount_u64(prefix) % 2 == 0) else -1


@njit(cache=True)
def decode_basis_jit(index, norbs, noccus, offsets, sizes):
    b = np.uint64(0)
    for n in range(len(norbs) - 1, -1, -1):
        d = sizes[n]
        sub_rank = index % d
        index //= d
        sub_bits = hash_decoder_jit(sub_rank, norbs[n], noccus[n])
        b |= np.uint64(sub_bits) << np.uint64(offsets[n])
    return b


@njit(cache=True)
def encode_basis_jit(b, norbs, noccus, offsets, strides, min_decode, max_decode):
    if b < min_decode or b > max_decode:
        return -1

    index = 0
    for n in range(len(norbs)):
        start = offsets[n]
        N = norbs[n]
        M = noccus[n]
        mask = (np.uint64(1) << np.uint64(N)) - np.uint64(1)
        sub_bits = (b >> np.uint64(start)) & mask

        if popcount_u64(sub_bits) != M:
            return -1

        index += hash_encoder_jit(sub_bits, N, M) * strides[n]

    return index


def decode_basis_py(index, norbs, noccus, offsets, sizes):
    norbs = np.asarray(norbs, dtype=np.int64)
    noccus = np.asarray(noccus, dtype=np.int64)
    offsets = np.asarray(offsets, dtype=np.int64)
    sizes = np.asarray(sizes, dtype=np.int64)
    return int(decode_basis_jit(np.int64(index), norbs, noccus, offsets, sizes))


def encode_basis_py(b, norbs, noccus, offsets, strides, min_decode, max_decode):
    norbs = np.asarray(norbs, dtype=np.int64)
    noccus = np.asarray(noccus, dtype=np.int64)
    offsets = np.asarray(offsets, dtype=np.int64)
    strides = np.asarray(strides, dtype=np.int64)
    return int(encode_basis_jit(np.uint64(b), norbs, noccus, offsets, strides,
                                np.uint64(min_decode), np.uint64(max_decode)))


@njit(cache=True)
def build_H_entries_lr(
    rb_norbs, rb_noccus, rb_offsets, rb_sizes, rb_strides, rb_min_decode, rb_max_decode,
    lb_norbs, lb_noccus, lb_offsets, lb_sizes, lb_strides, lb_min_decode, lb_max_decode,
    e_terms, e_vals, u_terms, u_vals, cstart, cend,
):
    """Build row/col/value triplets for an operator mapping rb -> lb.

    The right basis is used for decoding column configurations.
    The left basis is used for encoding the transformed configuration
    back into a row index.
    """
    rows = []
    cols = []
    vals = []

    norb_total = 0
    for x in rb_norbs:
        norb_total += x

    # 1-body terms: sum_ij e_{ij} f_i^dagger f_j
    for icfg in range(cstart, cend):
        b0 = decode_basis_jit(icfg, rb_norbs, rb_noccus, rb_offsets, rb_sizes)

        for t in range(e_terms.shape[0]):
            iorb = int(e_terms[t, 0])
            jorb = int(e_terms[t, 1])
            e = e_vals[t]

            if bit_at_orb_u64(b0, jorb, norb_total) == 0:
                continue
            s1 = sign_count_u64(b0, jorb, norb_total)
            b = set_zero_at_orb_u64(b0, jorb, norb_total)

            if bit_at_orb_u64(b, iorb, norb_total) == 1:
                continue
            s2 = sign_count_u64(b, iorb, norb_total)
            b = set_one_at_orb_u64(b, iorb, norb_total)

            jcfg = encode_basis_jit(
                b,
                lb_norbs, lb_noccus, lb_offsets, lb_strides, lb_min_decode, lb_max_decode,
            )
            if jcfg != -1:
                rows.append(jcfg)
                cols.append(icfg)
                vals.append(e * s1 * s2)

        # 2-body terms: sum_lkji u_{lkji} f_l^dagger f_k^dagger f_j f_i
        for t in range(u_terms.shape[0]):
            lorb = int(u_terms[t, 0])
            korb = int(u_terms[t, 1])
            jorb = int(u_terms[t, 2])
            iorb = int(u_terms[t, 3])
            u = u_vals[t]

            if iorb == jorb or korb == lorb:
                continue
            if bit_at_orb_u64(b0, iorb, norb_total) == 0:
                continue
            s1 = sign_count_u64(b0, iorb, norb_total)
            b = set_zero_at_orb_u64(b0, iorb, norb_total)

            if bit_at_orb_u64(b, jorb, norb_total) == 0:
                continue
            s2 = sign_count_u64(b, jorb, norb_total)
            b = set_zero_at_orb_u64(b, jorb, norb_total)

            if bit_at_orb_u64(b, korb, norb_total) == 1:
                continue
            s3 = sign_count_u64(b, korb, norb_total)
            b = set_one_at_orb_u64(b, korb, norb_total)

            if bit_at_orb_u64(b, lorb, norb_total) == 1:
                continue
            s4 = sign_count_u64(b, lorb, norb_total)
            b = set_one_at_orb_u64(b, lorb, norb_total)

            jcfg = encode_basis_jit(
                b,
                lb_norbs, lb_noccus, lb_offsets, lb_strides, lb_min_decode, lb_max_decode,
            )
            if jcfg != -1:
                rows.append(jcfg)
                cols.append(icfg)
                vals.append(u * s1 * s2 * s3 * s4)

    return (
        np.asarray(rows, dtype=np.int32),
        np.asarray(cols, dtype=np.int32),
        np.asarray(vals, dtype=np.complex128),
    )


def assemble_petsc_from_entries(comm, nl, nr, rows, cols, vals, nnz_guess_per_row=16):
    H = PETSc.Mat().create(comm=comm)
    H.setSizes(((None, nl), (None, nr)))
    H.setType(PETSc.Mat.Type.AIJ)
    # H.setType("aijcusparse")
    # H.setVecType(PETSc.Vec.Type.CUDA)
    H.setPreallocationNNZ(nnz_guess_per_row)
    H.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    H.setUp()

    print(f"H type in has script: {H.getType()}", flush=True)

    if rows.size:
        order = np.argsort(rows, kind="mergesort")
        rows = rows[order]
        cols = cols[order]
        vals = vals[order]

        start = 0
        while start < rows.size:
            r = int(rows[start])
            end = start + 1
            while end < rows.size and rows[end] == rows[start]:
                end += 1
            H.setValues(r, cols[start:end], vals[start:end], addv=PETSc.InsertMode.ADD_VALUES)
            start = end

    H.assemblyBegin()
    H.assemblyEnd()
    H_gpu = H.convert("aijcusparse")
    return H_gpu
    # return H


def get_H(comm, emat, umat, lb, rb=None, tol_e=1e-10, tol_u=1e-10, nnz_guess_per_row=None):
    """Assemble H = sum_ij emat_ij f_i^dagger f_j + sum_lkji umat_lkji f_l^dagger f_k^dagger f_j f_i.

    For rectangular operators, rb is the column basis and lb is the row basis.
    For square operators, omit rb and lb is used on both sides.
    """
    if rb is None:
        rb = lb

    nl, nr = lb.dim, rb.dim
    H = PETSc.Mat().create(comm=comm)
    H.setSizes(((None, nl), (None, nr)))
    H.setType(PETSc.Mat.Type.AIJ)

    if nnz_guess_per_row is None:
        ne = int(np.count_nonzero(np.abs(emat) > tol_e)) if emat is not None else 0
        nu = int(np.count_nonzero(np.abs(umat) > tol_u)) if umat is not None else 0
        nnz_guess_per_row = max(8, min(nr, 16 + ne + 2 * min(nu, 32)))

    H.setPreallocationNNZ(nnz_guess_per_row)
    H.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    H.setUp()

    cstart, cend = H.getOwnershipRangeColumn()
    print(f"cstart = {cstart}, cend = {cend}")

    rb_meta = rb.jit_args()
    lb_meta = lb.jit_args()

    if emat is not None:
        a1, a2 = np.nonzero(np.abs(emat) > tol_e)
        e_terms = np.stack((a1, a2), axis=-1).astype(np.int64)
        e_vals = emat[a1, a2].astype(np.complex128)
    else:
        e_terms = np.empty((0, 2), dtype=np.int64)
        e_vals = np.empty((0,), dtype=np.complex128)

    if umat is not None:
        a1, a2, a3, a4 = np.nonzero(np.abs(umat) > tol_u)
        u_terms = np.stack((a1, a2, a3, a4), axis=-1).astype(np.int64)
        u_vals = umat[a1, a2, a3, a4].astype(np.complex128)
    else:
        u_terms = np.empty((0, 4), dtype=np.int64)
        u_vals = np.empty((0,), dtype=np.complex128)

    rows, cols, vals = build_H_entries_lr(
        rb_meta[0], rb_meta[1], rb_meta[2], rb_meta[3], rb_meta[4], rb_meta[5], rb_meta[6],
        lb_meta[0], lb_meta[1], lb_meta[2], lb_meta[3], lb_meta[4], lb_meta[5], lb_meta[6],
        e_terms, e_vals, u_terms, u_vals, cstart, cend,
    )

    return assemble_petsc_from_entries(comm, nl, nr, rows, cols, vals, nnz_guess_per_row)
