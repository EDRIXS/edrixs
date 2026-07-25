"""Shared preparation and assertions for solver consistency checks."""

import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import ops


def exact_1v1c_reference_data(problem):
    """Diagonalize the tiny dense model used as a final-spectrum reference."""
    hmat_i, hmat_n, transitions = ops(*problem, backend="dense")
    eval_i, evec_i = np.linalg.eigh(hmat_i)
    eval_n, evec_n = np.linalg.eigh(hmat_n)
    transitions_eigenbasis = np.stack(
        [
            evec_n.conj().T @ transition @ evec_i
            for transition in transitions
        ]
    )
    return (
        hmat_i,
        hmat_n,
        transitions,
        eval_i,
        evec_i,
        eval_n,
        transitions_eigenbasis,
    )


def assert_problem_sparse_dense_equivalent(dense, sparse):
    """Compare backend-neutral setup outputs in dense and sparse-U forms."""
    assert_allclose(dense[0], sparse[0])
    assert_allclose(dense[3], sparse[3])
    assert dense[2].basis_int == sparse[2].basis_int
    assert dense[5].basis_int == sparse[5].basis_int
    assert_allclose(dense[6], sparse[6])

    for dense_u, sparse_u in ((dense[1], sparse[1]), (dense[4], sparse[4])):
        norbs = dense_u.shape[0]
        assert sp.isspmatrix_csr(sparse_u)
        assert sparse_u.shape == (norbs * norbs, norbs * norbs)
        assert_allclose(
            sparse_u.toarray(),
            dense_u.reshape(norbs * norbs, norbs * norbs),
            atol=1e-13,
        )


def assert_dense_and_scipy_ops_match(dense_problem, sparse_problem, seed=0):
    """Compare dense and SciPy operator actions produced from the same model."""
    hmat_i, hmat_n, transitions = ops(*dense_problem, backend="dense")
    hmat_i_sp, hmat_n_sp, transitions_sp = ops(
        *sparse_problem,
        backend="scipy",
    )

    rng = np.random.default_rng(seed)
    vec_i = rng.normal(size=hmat_i.shape[1]) + 1j * rng.normal(
        size=hmat_i.shape[1]
    )
    vec_n = rng.normal(size=hmat_n.shape[1]) + 1j * rng.normal(
        size=hmat_n.shape[1]
    )

    assert_allclose(hmat_i_sp @ vec_i, hmat_i @ vec_i, atol=2e-12)
    assert_allclose(hmat_i_sp.H @ vec_i, hmat_i.conj().T @ vec_i, atol=2e-12)
    assert_allclose(hmat_n_sp @ vec_n, hmat_n @ vec_n, atol=2e-12)
    assert_allclose(hmat_n_sp.H @ vec_n, hmat_n.conj().T @ vec_n, atol=2e-12)

    for transition_sp, transition in zip(transitions_sp, transitions):
        assert_allclose(transition_sp @ vec_i, transition @ vec_i, atol=2e-12)
        assert_allclose(
            transition_sp.H @ vec_n,
            transition.conj().T @ vec_n,
            atol=2e-12,
        )


def assert_only_transition_block_is_populated(trans_mat, row_slice, col_slice):
    """Check that setup places photon transitions only in the selected shell."""
    expected = np.zeros_like(trans_mat)
    expected[:, row_slice, col_slice] = trans_mat[:, row_slice, col_slice]
    assert_allclose(trans_mat, expected)
    assert np.linalg.norm(trans_mat[:, row_slice, col_slice]) > 0


def siam_kwargs(siam_type):
    """Return tiny equivalent SIAM inputs for both supported representations."""
    common = dict(
        shell_name=("s", "p"),
        nbath=1,
        siam_type=siam_type,
        v_noccu=1,
        static_core_pot=0.35,
        c_level=-2.7,
        c_soc=0.06,
        slater=([0.3], [0.45]),
    )

    if siam_type == 0:
        common.update(
            imp_mat=np.array([[0.12, 0.02j], [-0.02j, -0.12]]),
            imp_mat_n=np.array([[0.18, -0.01j], [0.01j, -0.08]]),
            bath_level=np.array([[0.4, 0.55]]),
            bath_level_n=np.array([[0.45, 0.5]]),
            hyb=np.array([[0.16, 0.08j]]),
            hyb_n=np.array([[0.19, -0.04j]]),
        )
    else:
        hopping = np.array(
            [
                [0.10, 0.02j, 0.15, 0.00],
                [-0.02j, -0.10, 0.00, -0.12],
                [0.15, 0.00, 0.45, 0.01j],
                [0.00, -0.12, -0.01j, 0.52],
            ],
            dtype=complex,
        )
        hopping_n = hopping.copy()
        hopping_n[0, 0] += 0.07
        hopping_n[1, 1] -= 0.03
        common.update(hopping=hopping, hopping_n=hopping_n)

    return common


def random_hermitian(rng, size):
    """Build a tiny Hermitian matrix for backend-independent spectrum checks."""
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return (raw + raw.conj().T) / 2
