import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import ops, setup_2v1c, setup_siam


def _assert_problem_sparse_dense_equivalent(dense, sparse):
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


def _assert_ops_match(dense_problem, sparse_problem, seed=0):
    h_i, h_n, transitions = ops(*dense_problem, backend="dense")
    h_i_sp, h_n_sp, transitions_sp = ops(*sparse_problem, backend="scipy")

    rng = np.random.default_rng(seed)
    vec_i = rng.normal(size=h_i.shape[1]) + 1j * rng.normal(size=h_i.shape[1])
    vec_n = rng.normal(size=h_n.shape[1]) + 1j * rng.normal(size=h_n.shape[1])

    assert_allclose(h_i_sp @ vec_i, h_i @ vec_i, atol=2e-12)
    assert_allclose(h_i_sp.H @ vec_i, h_i.conj().T @ vec_i, atol=2e-12)
    assert_allclose(h_n_sp @ vec_n, h_n @ vec_n, atol=2e-12)
    assert_allclose(h_n_sp.H @ vec_n, h_n.conj().T @ vec_n, atol=2e-12)

    for transition_sp, transition in zip(transitions_sp, transitions):
        assert_allclose(transition_sp @ vec_i, transition @ vec_i, atol=2e-12)
        assert_allclose(
            transition_sp.H @ vec_n,
            transition.conj().T @ vec_n,
            atol=2e-12,
        )


def _assert_only_transition_block_is_populated(trans_mat, row_slice, col_slice):
    expected = np.zeros_like(trans_mat)
    expected[:, row_slice, col_slice] = trans_mat[:, row_slice, col_slice]
    assert_allclose(trans_mat, expected)
    assert np.linalg.norm(trans_mat[:, row_slice, col_slice]) > 0


@pytest.mark.integration
@pytest.mark.parametrize("trans_to_which", [1, 2])
def test_setup_2v1c_sparse_dense_and_ops_equivalence(trans_to_which):
    # Both valence shells are s-like and the p core can transition into either.
    hopping = np.array(
        [
            [0.04 + 0.02j, -0.07j],
            [0.03, -0.02 + 0.01j],
        ],
        dtype=complex,
    )
    kwargs = dict(
        shell_name=("s", "s", "p"),
        shell_level=(0.15, -0.25, -3.5),
        c_soc=0.08,
        v_tot_noccu=1,
        slater=([0.4], [0.55]),
        v1_othermat=np.array([[0.03, 0.01], [0.01, -0.03]]),
        v2_othermat=np.array([[-0.02, 0.015j], [-0.015j, 0.02]]),
        hopping_v1v2=hopping,
        trans_to_which=trans_to_which,
    )

    dense = setup_2v1c(**kwargs, sparse_U=False)
    sparse = setup_2v1c(**kwargs, sparse_U=True)

    _assert_problem_sparse_dense_equivalent(dense, sparse)
    _assert_ops_match(dense, sparse, seed=10 + trans_to_which)

    # For (s, s, p), v1 and v2 each contain two spin-orbitals and
    # the core begins at orbital four.
    if trans_to_which == 1:
        _assert_only_transition_block_is_populated(dense[6], slice(0, 2), slice(4, 10))
    else:
        _assert_only_transition_block_is_populated(dense[6], slice(2, 4), slice(4, 10))

    assert_allclose(dense[0][2:4, 0:2], hopping.conj().T)
    assert_allclose(dense[0][0:2, 2:4], hopping)
    assert_allclose(dense[0], dense[0].conj().T)
    assert_allclose(dense[3], dense[3].conj().T)


def test_setup_2v1c_rejects_unknown_transition_target():
    with pytest.raises(Exception, match="Unknown trans_to_which"):
        setup_2v1c(("s", "s", "p"), trans_to_which=3)


def _siam_kwargs(siam_type):
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


@pytest.mark.integration
@pytest.mark.parametrize("siam_type", [0, 1])
def test_setup_siam_sparse_dense_and_ops_equivalence(siam_type):
    kwargs = _siam_kwargs(siam_type)
    dense = setup_siam(**kwargs, sparse_U=False)
    sparse = setup_siam(**kwargs, sparse_U=True)

    _assert_problem_sparse_dense_equivalent(dense, sparse)
    _assert_ops_match(dense, sparse, seed=20 + siam_type)

    # For (s, p) with one bath, orbitals 0:2 are the impurity, 2:4 the bath,
    # and 4:10 the core. The photon transition must not act on the bath.
    transition = dense[6]
    _assert_only_transition_block_is_populated(transition, slice(0, 2), slice(4, 10))
    assert_allclose(transition[:, 2:4, :], 0.0)


@pytest.mark.integration
def test_setup_siam_static_core_potential_only_shifts_intermediate_impurity():
    base = setup_siam(
        ("s", "p"),
        1,
        v_noccu=1,
        static_core_pot=0.0,
        c_level=0.0,
        sparse_U=False,
    )
    shifted = setup_siam(
        ("s", "p"),
        1,
        v_noccu=1,
        static_core_pot=0.6,
        c_level=0.0,
        sparse_U=False,
    )

    assert_allclose(shifted[0], base[0])
    expected = np.zeros_like(base[3])
    expected[0:2, 0:2] = -0.6 * np.eye(2)
    assert_allclose(shifted[3] - base[3], expected)


def test_setup_siam_rejects_unknown_type():
    with pytest.raises(Exception, match="Unknown siam_type"):
        setup_siam(("s", "p"), 1, siam_type=7)
