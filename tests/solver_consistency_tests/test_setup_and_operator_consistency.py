"""Consistency checks for model setup and many-body operator construction.

These checks follow the command chain from ``setup_*`` through ``get_ops`` but stop
before any eigensolver or spectral solver.  They compare alternative dense and
SciPy sparse representations of the same physical model.
"""

import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.models import setup_1v1c, setup_2v1c, setup_siam
from edrixs.solvers import get_ops

from _helpers import (
    assert_dense_and_scipy_get_ops_match,
    assert_only_transition_block_is_populated,
    assert_problem_sparse_dense_equivalent,
    siam_kwargs,
)

pytestmark = pytest.mark.integration


@pytest.mark.parametrize("trans_to_which", [1, 2])
def test_setup_2v1c_sparse_dense_and_get_ops_equivalence(trans_to_which):
    """Compare both representations through the 2v1c setup-to-``get_ops`` chain.

    The check builds one model with dense and flattened sparse interactions,
    verifies the same orbital data and transition target, and confirms that
    dense and SciPy many-body operators act identically before spectroscopy.
    """
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

    assert_problem_sparse_dense_equivalent(dense, sparse)
    assert_dense_and_scipy_get_ops_match(dense, sparse, seed=10 + trans_to_which)

    if trans_to_which == 1:
        assert_only_transition_block_is_populated(
            dense[6],
            slice(0, 2),
            slice(4, 10),
        )
    else:
        assert_only_transition_block_is_populated(
            dense[6],
            slice(2, 4),
            slice(4, 10),
        )

    assert_allclose(dense[0][2:4, 0:2], hopping.conj().T)
    assert_allclose(dense[0][0:2, 2:4], hopping)
    assert_allclose(dense[0], dense[0].conj().T)
    assert_allclose(dense[3], dense[3].conj().T)


@pytest.mark.parametrize("siam_type", [0, 1])
def test_setup_siam_sparse_dense_and_get_ops_equivalence(siam_type):
    """Compare both representations through the SIAM setup-to-``get_ops`` chain.

    The check covers both supported SIAM input forms, confirms identical dense
    and sparse model data and operator actions, and verifies that the photon
    transition reaches the impurity but not bath orbitals.
    """
    kwargs = siam_kwargs(siam_type)
    dense = setup_siam(**kwargs, sparse_U=False)
    sparse = setup_siam(**kwargs, sparse_U=True)

    assert_problem_sparse_dense_equivalent(dense, sparse)
    assert_dense_and_scipy_get_ops_match(dense, sparse, seed=20 + siam_type)

    transition = dense[6]
    assert_only_transition_block_is_populated(
        transition,
        slice(0, 2),
        slice(4, 10),
    )
    assert_allclose(transition[:, 2:4, :], 0.0)


def test_setup_siam_static_core_potential_only_shifts_intermediate_impurity():
    """Check the SIAM core-hole term before operator construction and spectra.

    The static core potential should leave the initial orbital Hamiltonian
    unchanged and shift only the impurity block of the intermediate model that
    later enters XAS and the RIXS correction-vector solve.
    """
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


def test_setup_1v1c_sparse_u_matches_dense_u(small_1v1c_kwargs):
    """Compare dense and sparse interaction outputs from ``setup_1v1c``.

    This check verifies the first stage of the 1v1c SciPy path: model definition
    must be unchanged when Coulomb tensors are stored sparsely for the later
    ``get_ops(..., backend="scipy")`` construction.
    """
    dense = setup_1v1c(**small_1v1c_kwargs, sparse_U=False)
    sparse = setup_1v1c(**small_1v1c_kwargs, sparse_U=True)

    for dense_u, sparse_u in ((dense[1], sparse[1]), (dense[4], sparse[4])):
        norbs = dense_u.shape[0]
        assert sp.isspmatrix_csr(sparse_u)
        assert_allclose(
            sparse_u.toarray(),
            dense_u.reshape(norbs * norbs, norbs * norbs),
        )

    assert dense[2].basis_int == sparse[2].basis_int
    assert dense[5].basis_int == sparse[5].basis_int
    assert_allclose(dense[0], sparse[0])
    assert_allclose(dense[3], sparse[3])
    assert_allclose(dense[6], sparse[6])


def test_ops_scipy_backend_matches_dense_operator_action(small_1v1c_problem):
    """Compare dense and SciPy operators at the setup-to-solver boundary.

    Starting from one backend-neutral 1v1c model, this check confirms that the
    two ``get_ops`` backends produce Hamiltonian and transition actions that agree
    before either representation is used by ED, XAS, or RIXS.
    """
    hmat_i_sp, hmat_n_sp, trans_sp = get_ops(
        *small_1v1c_problem,
        backend="scipy",
    )
    hmat_i, hmat_n, transitions = get_ops(
        *small_1v1c_problem,
        backend="dense",
    )
    rng = np.random.default_rng(5)

    vec_i = rng.normal(size=hmat_i.shape[1]) + 1j * rng.normal(
        size=hmat_i.shape[1]
    )
    vec_n = rng.normal(size=hmat_n.shape[1]) + 1j * rng.normal(
        size=hmat_n.shape[1]
    )
    assert_allclose(hmat_i_sp @ vec_i, hmat_i @ vec_i, atol=1e-12)
    assert_allclose(
        hmat_i_sp.conj().T @ vec_i,
        hmat_i.conj().T @ vec_i,
        atol=1e-12,
    )
    assert_allclose(hmat_n_sp @ vec_n, hmat_n @ vec_n, atol=1e-12)

    for sparse_op, dense_op in zip(trans_sp, transitions):
        assert_allclose(sparse_op @ vec_i, dense_op @ vec_i, atol=1e-12)
        assert_allclose(
            sparse_op.conj().T @ vec_n,
            dense_op.conj().T @ vec_n,
            atol=1e-12,
        )
