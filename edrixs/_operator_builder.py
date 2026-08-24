"""Backend-independent kernels for many-body operator matrix entries."""

from __future__ import annotations

import numpy as np

from .fock_basis import FockBasis, FockBinByN


_NUMBA_KERNELS = None


def _operator_terms(emat, umat, norbs, tol_e, tol_u):
    if emat is not None:
        emat = np.asarray(emat)
        if emat.shape != (norbs, norbs):
            raise ValueError(
                "emat has shape {}, expected {}".format(
                    emat.shape, (norbs, norbs)
                )
            )
        iorb, jorb = np.nonzero(np.abs(emat) > tol_e)
        e_terms = np.stack((iorb, jorb), axis=-1).astype(np.int64)
        e_vals = emat[iorb, jorb].astype(np.complex128)
    else:
        e_terms = np.empty((0, 2), dtype=np.int64)
        e_vals = np.empty((0,), dtype=np.complex128)

    if umat is None:
        u_terms = np.empty((0, 4), dtype=np.int64)
        u_vals = np.empty((0,), dtype=np.complex128)
    elif hasattr(umat, 'tocoo') and len(umat.shape) == 2:
        expected = (norbs * norbs, norbs * norbs)
        if umat.shape != expected:
            raise ValueError(
                "sparse umat has shape {}, expected {}".format(umat.shape, expected)
            )
        coo = umat.tocoo()
        keep = np.abs(coo.data) > tol_u
        rows = np.asarray(coo.row[keep], dtype=np.int64)
        cols = np.asarray(coo.col[keep], dtype=np.int64)
        lorb, korb = np.divmod(rows, norbs)
        jorb, iorb = np.divmod(cols, norbs)
        u_terms = np.stack((lorb, korb, jorb, iorb), axis=-1)
        u_vals = np.asarray(coo.data[keep], dtype=np.complex128)
    else:
        umat = np.asarray(umat)
        expected = (norbs, norbs, norbs, norbs)
        if umat.shape != expected:
            raise ValueError(
                "dense umat has shape {}, expected {}".format(umat.shape, expected)
            )
        lorb, korb, jorb, iorb = np.nonzero(np.abs(umat) > tol_u)
        u_terms = np.stack((lorb, korb, jorb, iorb), axis=-1).astype(np.int64)
        u_vals = umat[lorb, korb, jorb, iorb].astype(np.complex128)

    return e_terms, e_vals, u_terms, u_vals


def _sign_count(state, orbital, norbs):
    prefix = int(state) >> (norbs - orbital)
    return 1 if prefix.bit_count() % 2 == 0 else -1


def _build_entries_python(lb, rb, e_terms, e_vals, u_terms, u_vals, cstart, cend):
    rows = []
    cols = []
    vals = []
    norbs = rb.norbs

    for column in range(cstart, cend):
        state0 = rb.decode(column)

        for term, value in zip(e_terms, e_vals):
            iorb, jorb = int(term[0]), int(term[1])
            bit_j = 1 << (norbs - 1 - jorb)
            if not state0 & bit_j:
                continue
            sign_1 = _sign_count(state0, jorb, norbs)
            state = state0 ^ bit_j

            bit_i = 1 << (norbs - 1 - iorb)
            if state & bit_i:
                continue
            sign_2 = _sign_count(state, iorb, norbs)
            state |= bit_i

            try:
                row = lb.encode(state)
            except KeyError:
                continue
            rows.append(row)
            cols.append(column)
            vals.append(value * sign_1 * sign_2)

        for term, value in zip(u_terms, u_vals):
            lorb, korb, jorb, iorb = map(int, term)
            if iorb == jorb or korb == lorb:
                continue

            bit_i = 1 << (norbs - 1 - iorb)
            if not state0 & bit_i:
                continue
            sign_1 = _sign_count(state0, iorb, norbs)
            state = state0 ^ bit_i

            bit_j = 1 << (norbs - 1 - jorb)
            if not state & bit_j:
                continue
            sign_2 = _sign_count(state, jorb, norbs)
            state ^= bit_j

            bit_k = 1 << (norbs - 1 - korb)
            if state & bit_k:
                continue
            sign_3 = _sign_count(state, korb, norbs)
            state |= bit_k

            bit_l = 1 << (norbs - 1 - lorb)
            if state & bit_l:
                continue
            sign_4 = _sign_count(state, lorb, norbs)
            state |= bit_l

            try:
                row = lb.encode(state)
            except KeyError:
                continue
            rows.append(row)
            cols.append(column)
            vals.append(value * sign_1 * sign_2 * sign_3 * sign_4)

    return (
        np.asarray(rows, dtype=np.int64),
        np.asarray(cols, dtype=np.int64),
        np.asarray(vals, dtype=np.complex128),
    )


def _get_numba_kernels():
    global _NUMBA_KERNELS
    if _NUMBA_KERNELS is not None:
        return _NUMBA_KERNELS

    try:
        from numba import njit
    except ImportError as exc:
        raise ImportError(
            "use_numba=True requires numba; install numba or disable JIT construction"
        ) from exc

    @njit(inline='always')
    def comb_jit(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n - k:
            k = n - k
        result = 1
        for i in range(k):
            result = (result * (n - i)) // (i + 1)
        return result

    @njit(inline='always')
    def popcount_u64(value):
        count = 0
        while value != np.uint64(0):
            value &= value - np.uint64(1)
            count += 1
        return count

    @njit(inline='always')
    def hash_decoder_jit(rank, norb, nocc):
        state = np.uint64(0)
        j = nocc
        i = norb - 1
        c = comb_jit(i, j)
        while i >= 0 and j > 0:
            if c <= rank:
                state |= np.uint64(1) << np.uint64(i)
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

    @njit(inline='always')
    def hash_encoder_jit(state, norb):
        if norb < 64:
            state &= (np.uint64(1) << np.uint64(norb)) - np.uint64(1)
        rank = 0
        k = 1
        while state != np.uint64(0):
            lsb = state & (np.uint64(0) - state)
            pos = 0
            tmp = lsb
            while tmp > np.uint64(1):
                tmp >>= np.uint64(1)
                pos += 1
            rank += comb_jit(pos, k)
            k += 1
            state ^= lsb
        return rank

    @njit(inline='always')
    def sign_count_u64(state, orbital, norbs):
        bitpos = np.uint64(norbs - 1 - orbital)
        if bitpos == np.uint64(63):
            prefix = np.uint64(0)
        else:
            prefix = state >> (bitpos + np.uint64(1))
        return 1 if popcount_u64(prefix) % 2 == 0 else -1

    @njit(inline='always')
    def decode_combinadic_jit(index, norbs, noccus, offsets, sizes):
        state = np.uint64(0)
        for n in range(len(norbs)):
            size = sizes[n]
            shell_rank = index % size
            index //= size
            colex_rank = size - 1 - shell_rank
            shell_state = hash_decoder_jit(colex_rank, norbs[n], noccus[n])
            state |= shell_state << np.uint64(offsets[n])
        return state

    @njit(inline='always')
    def encode_combinadic_jit(state, norbs, noccus, offsets, sizes, strides):
        index = 0
        for n in range(len(norbs)):
            norb = norbs[n]
            if norb == 64:
                mask = np.uint64(0xffffffffffffffff)
            else:
                mask = (np.uint64(1) << np.uint64(norb)) - np.uint64(1)
            shell_state = (state >> np.uint64(offsets[n])) & mask
            if popcount_u64(shell_state) != noccus[n]:
                return -1
            colex_rank = hash_encoder_jit(shell_state, norb)
            shell_rank = sizes[n] - 1 - colex_rank
            index += shell_rank * strides[n]
        return index

    @njit
    def build_combinadic(
        rb_norbs, rb_noccus, rb_offsets, rb_sizes, rb_strides,
        lb_norbs, lb_noccus, lb_offsets, lb_sizes, lb_strides,
        e_terms, e_vals, u_terms, u_vals, cstart, cend,
    ):
        rows = []
        cols = []
        vals = []
        norbs = 0
        for value in rb_norbs:
            norbs += value

        for column in range(cstart, cend):
            state0 = decode_combinadic_jit(
                column, rb_norbs, rb_noccus, rb_offsets, rb_sizes
            )

            for t in range(e_terms.shape[0]):
                iorb = int(e_terms[t, 0])
                jorb = int(e_terms[t, 1])
                bit_j = np.uint64(1) << np.uint64(norbs - 1 - jorb)
                if state0 & bit_j == 0:
                    continue
                sign_1 = sign_count_u64(state0, jorb, norbs)
                state = state0 ^ bit_j

                bit_i = np.uint64(1) << np.uint64(norbs - 1 - iorb)
                if state & bit_i != 0:
                    continue
                sign_2 = sign_count_u64(state, iorb, norbs)
                state |= bit_i

                row = encode_combinadic_jit(
                    state, lb_norbs, lb_noccus, lb_offsets, lb_sizes, lb_strides
                )
                if row != -1:
                    rows.append(row)
                    cols.append(column)
                    vals.append(e_vals[t] * sign_1 * sign_2)

            for t in range(u_terms.shape[0]):
                lorb = int(u_terms[t, 0])
                korb = int(u_terms[t, 1])
                jorb = int(u_terms[t, 2])
                iorb = int(u_terms[t, 3])
                if iorb == jorb or korb == lorb:
                    continue

                bit_i = np.uint64(1) << np.uint64(norbs - 1 - iorb)
                if state0 & bit_i == 0:
                    continue
                sign_1 = sign_count_u64(state0, iorb, norbs)
                state = state0 ^ bit_i

                bit_j = np.uint64(1) << np.uint64(norbs - 1 - jorb)
                if state & bit_j == 0:
                    continue
                sign_2 = sign_count_u64(state, jorb, norbs)
                state ^= bit_j

                bit_k = np.uint64(1) << np.uint64(norbs - 1 - korb)
                if state & bit_k != 0:
                    continue
                sign_3 = sign_count_u64(state, korb, norbs)
                state |= bit_k

                bit_l = np.uint64(1) << np.uint64(norbs - 1 - lorb)
                if state & bit_l != 0:
                    continue
                sign_4 = sign_count_u64(state, lorb, norbs)
                state |= bit_l

                row = encode_combinadic_jit(
                    state, lb_norbs, lb_noccus, lb_offsets, lb_sizes, lb_strides
                )
                if row != -1:
                    rows.append(row)
                    cols.append(column)
                    vals.append(u_vals[t] * sign_1 * sign_2 * sign_3 * sign_4)

        return (
            np.asarray(rows, dtype=np.int64),
            np.asarray(cols, dtype=np.int64),
            np.asarray(vals, dtype=np.complex128),
        )

    @njit
    def build_explicit(
        rb_states, lb_lookup, norbs,
        e_terms, e_vals, u_terms, u_vals, cstart, cend,
    ):
        rows = []
        cols = []
        vals = []

        for column in range(cstart, cend):
            state0 = rb_states[column]

            for t in range(e_terms.shape[0]):
                iorb = int(e_terms[t, 0])
                jorb = int(e_terms[t, 1])
                bit_j = np.uint64(1) << np.uint64(norbs - 1 - jorb)
                if state0 & bit_j == 0:
                    continue
                sign_1 = sign_count_u64(state0, jorb, norbs)
                state = state0 ^ bit_j

                bit_i = np.uint64(1) << np.uint64(norbs - 1 - iorb)
                if state & bit_i != 0:
                    continue
                sign_2 = sign_count_u64(state, iorb, norbs)
                state |= bit_i

                if state in lb_lookup:
                    row = lb_lookup[state]
                    rows.append(row)
                    cols.append(column)
                    vals.append(e_vals[t] * sign_1 * sign_2)

            for t in range(u_terms.shape[0]):
                lorb = int(u_terms[t, 0])
                korb = int(u_terms[t, 1])
                jorb = int(u_terms[t, 2])
                iorb = int(u_terms[t, 3])
                if iorb == jorb or korb == lorb:
                    continue

                bit_i = np.uint64(1) << np.uint64(norbs - 1 - iorb)
                if state0 & bit_i == 0:
                    continue
                sign_1 = sign_count_u64(state0, iorb, norbs)
                state = state0 ^ bit_i

                bit_j = np.uint64(1) << np.uint64(norbs - 1 - jorb)
                if state & bit_j == 0:
                    continue
                sign_2 = sign_count_u64(state, jorb, norbs)
                state ^= bit_j

                bit_k = np.uint64(1) << np.uint64(norbs - 1 - korb)
                if state & bit_k != 0:
                    continue
                sign_3 = sign_count_u64(state, korb, norbs)
                state |= bit_k

                bit_l = np.uint64(1) << np.uint64(norbs - 1 - lorb)
                if state & bit_l != 0:
                    continue
                sign_4 = sign_count_u64(state, lorb, norbs)
                state |= bit_l

                if state in lb_lookup:
                    row = lb_lookup[state]
                    rows.append(row)
                    cols.append(column)
                    vals.append(u_vals[t] * sign_1 * sign_2 * sign_3 * sign_4)

        return (
            np.asarray(rows, dtype=np.int64),
            np.asarray(cols, dtype=np.int64),
            np.asarray(vals, dtype=np.complex128),
        )

    _NUMBA_KERNELS = build_combinadic, build_explicit
    return _NUMBA_KERNELS


def _prepare_entries_kernel(
        lb, rb, e_terms, e_vals, u_terms, u_vals, *, use_numba=False):
    """Prepare a reusable entry kernel for backend-selected column ranges."""
    if not use_numba:
        def build_range(cstart, cend):
            return _build_entries_python(
                lb, rb, e_terms, e_vals, u_terms, u_vals, cstart, cend
            )

        return build_range

    if lb.norbs > 64 or rb.norbs > 64:
        raise ValueError("Numba operator construction supports at most 64 orbitals")

    build_combinadic, build_explicit = _get_numba_kernels()

    if isinstance(lb, FockBinByN) and isinstance(rb, FockBinByN):
        rb_meta = rb.jit_args()
        lb_meta = lb.jit_args()

        def build_range(cstart, cend):
            return build_combinadic(
                rb_meta[0], rb_meta[1], rb_meta[2], rb_meta[3], rb_meta[4],
                lb_meta[0], lb_meta[1], lb_meta[2], lb_meta[3], lb_meta[4],
                e_terms, e_vals, u_terms, u_vals, cstart, cend,
            )

        return build_range

    if isinstance(lb, FockBasis) and isinstance(rb, FockBasis):
        try:
            from numba import types
            from numba.typed import Dict
        except ImportError as exc:
            raise ImportError(
                "use_numba=True requires numba; install numba or disable JIT construction"
            ) from exc

        rb_states = np.asarray(rb.basis_int, dtype=np.uint64)
        lb_lookup = Dict.empty(key_type=types.uint64, value_type=types.int64)
        for index, state in enumerate(lb.basis_int):
            lb_lookup[np.uint64(state)] = np.int64(index)

        def build_range(cstart, cend):
            return build_explicit(
                rb_states, lb_lookup, rb.norbs,
                e_terms, e_vals, u_terms, u_vals, cstart, cend,
            )

        return build_range

    raise TypeError(
        "Numba construction requires matching explicit or combinadic basis representations"
    )


def prepare_operator_entry_kernel(
        emat, umat, lb, rb=None, *, tol_e=1e-10, tol_u=1e-10,
        use_numba=False):
    """Prepare a reusable kernel that generates entries for requested columns.

    The caller owns the work decomposition and matrix assembly policy. Calling
    the returned function as ``build_range(cstart, cend)`` produces only the
    contributions for that column interval.
    """
    if rb is None:
        rb = lb
    if lb.norbs != rb.norbs:
        raise ValueError("left and right Fock bases must have the same norbs")

    e_terms, e_vals, u_terms, u_vals = _operator_terms(
        emat, umat, rb.norbs, tol_e, tol_u
    )
    return _prepare_entries_kernel(
        lb, rb, e_terms, e_vals, u_terms, u_vals, use_numba=use_numba
    )


def build_operator_entries(
        emat, umat, lb, rb=None, *, tol_e=1e-10, tol_u=1e-10,
        cstart=0, cend=None, use_numba=False):
    """Return entries for one requested column range.

    This compatibility wrapper leaves the work-range choice with the caller;
    backends that need repeated bounded ranges should prepare the reusable
    builder with :func:`prepare_operator_entry_kernel` instead.
    """
    if rb is None:
        rb = lb
    if cend is None:
        cend = len(rb)

    build_range = prepare_operator_entry_kernel(
        emat, umat, lb, rb, tol_e=tol_e, tol_u=tol_u, use_numba=use_numba
    )
    return build_range(cstart, cend)
