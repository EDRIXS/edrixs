import multiprocessing as mp
import numpy as np
import pytest
import scipy.sparse as sp
from numpy.testing import assert_allclose

from edrixs.solvers import (
    ed_krylov_scipy,
    ops,
    rixs_1v1c_py,
    rixs_krylov_scipy,
    setup_1v1c,
    xas_1v1c_py,
    xas_krylov_scipy,
)


def _exact_problem_data(problem):
    h_i, h_n, trans = ops(*problem, backend="dense")
    eval_i, evec_i = np.linalg.eigh(h_i)
    eval_n, evec_n = np.linalg.eigh(h_n)
    trans_eigenbasis = np.stack(
        [evec_n.conj().T @ operator @ evec_i for operator in trans]
    )
    return h_i, h_n, trans, eval_i, evec_i, eval_n, trans_eigenbasis


@pytest.mark.integration
def test_setup_sparse_u_matches_dense_u(small_1v1c_kwargs):
    dense = setup_1v1c(**small_1v1c_kwargs, sparse_U=False)
    sparse = setup_1v1c(**small_1v1c_kwargs, sparse_U=True)

    for dense_u, sparse_u in ((dense[1], sparse[1]), (dense[4], sparse[4])):
        n = dense_u.shape[0]
        assert sp.isspmatrix_csr(sparse_u)
        assert_allclose(sparse_u.toarray(), dense_u.reshape(n * n, n * n))

    assert dense[2].basis_int == sparse[2].basis_int
    assert dense[5].basis_int == sparse[5].basis_int
    assert_allclose(dense[0], sparse[0])
    assert_allclose(dense[3], sparse[3])
    assert_allclose(dense[6], sparse[6])


@pytest.mark.integration
def test_ops_scipy_backend_matches_dense_operator_action(small_1v1c_problem):
    h_i_sp, h_n_sp, trans_sp = ops(*small_1v1c_problem, backend="scipy")
    h_i, h_n, trans = ops(*small_1v1c_problem, backend="dense")
    rng = np.random.default_rng(5)

    vi = rng.normal(size=h_i.shape[1]) + 1j * rng.normal(size=h_i.shape[1])
    vn = rng.normal(size=h_n.shape[1]) + 1j * rng.normal(size=h_n.shape[1])
    assert_allclose(h_i_sp @ vi, h_i @ vi, atol=1e-12)
    assert_allclose(h_i_sp.H @ vi, h_i.conj().T @ vi, atol=1e-12)
    assert_allclose(h_n_sp @ vn, h_n @ vn, atol=1e-12)

    for sparse_op, dense_op in zip(trans_sp, trans):
        assert_allclose(sparse_op @ vi, dense_op @ vi, atol=1e-12)
        assert_allclose(sparse_op.H @ vn, dense_op.conj().T @ vn, atol=1e-12)


def test_ed_krylov_returns_lowest_sorted_eigenpairs():
    diagonal = np.array([-3.0, -1.0, 0.5, 1.0, 2.0, 4.0, 5.0, 8.0])
    h = sp.diags(diagonal, format="csr")

    evals, evecs = ed_krylov_scipy(
        h,
        num_gs=2,
        blocksize=2,
        initial_guess=np.eye(len(diagonal), 2),
        tol=1e-12,
    )

    assert_allclose(evals, diagonal[:2], atol=1e-12)
    assert_allclose(h @ evecs, evecs * evals, atol=1e-11)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"num_gs": 0}, "positive"),
        ({"num_gs": 2, "blocksize": 1}, "greater than or equal"),
        ({"num_gs": 1, "blocksize": 9}, "cannot exceed"),
        (
            {"num_gs": 1, "blocksize": 1, "initial_guess": np.ones((7, 1))},
            "initial_guess",
        ),
    ],
)
def test_ed_krylov_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        ed_krylov_scipy(np.eye(8), **kwargs)


@pytest.mark.integration
def test_xas_krylov_matches_existing_dense_solver(small_1v1c_problem):
    _, h_n, trans, eval_i, evec_i, eval_n, trans_eig = _exact_problem_data(
        small_1v1c_problem
    )
    kept = [0]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.linspace(center - 0.8, center + 0.8, 9)
    pol_type = [("linear", 0.2), ("isotropic", 0.0)]

    sparse_result = xas_krylov_scipy(
        eval_i[kept],
        evec_i[:, kept],
        h_n,
        trans,
        ominc,
        gamma_c=0.25,
        thin=0.7,
        phi=0.1,
        pol_type=pol_type,
        temperature=20.0,
        nkryl=h_n.shape[0],
    )
    dense_result = xas_1v1c_py(
        eval_i,
        eval_n,
        trans_eig,
        ominc,
        gamma_c=0.25,
        thin=0.7,
        phi=0.1,
        pol_type=pol_type,
        gs_list=kept,
        temperature=20.0,
    )

    assert_allclose(sparse_result, dense_result, rtol=2e-8, atol=2e-10)


@pytest.mark.integration
def test_rixs_krylov_matches_existing_dense_solver(small_1v1c_problem):
    h_i, h_n, trans, eval_i, evec_i, eval_n, trans_eig = _exact_problem_data(
        small_1v1c_problem
    )
    kept = [0]
    center = float(np.median(eval_n) - eval_i[0])
    ominc = np.array([center - 0.2, center + 0.25])
    eloss = np.linspace(-0.4, 1.2, 11)
    pol_type = [("linear", 0.1, "linear", -0.3)]

    sparse_result = rixs_krylov_scipy(
        eval_i[kept],
        evec_i[:, kept],
        h_i,
        h_n,
        trans,
        ominc,
        eloss,
        gamma_c=0.25,
        gamma_f=0.08,
        thin=0.7,
        thout=1.1,
        phi=0.2,
        pol_type=pol_type,
        temperature=20.0,
        nkryl=h_i.shape[0],
        linsys_tol=1e-12,
        linsys_maxiter=500,
        linsys_restart=max(5, h_n.shape[0]),
        workers=1,
    )
    dense_result = rixs_1v1c_py(
        eval_i,
        eval_n,
        trans_eig,
        ominc,
        eloss,
        gamma_c=0.25,
        gamma_f=0.08,
        thin=0.7,
        thout=1.1,
        phi=0.2,
        pol_type=pol_type,
        gs_list=kept,
        temperature=20.0,
    )

    assert_allclose(sparse_result, dense_result, rtol=2e-7, atol=2e-9)


@pytest.mark.integration
def test_rixs_return_poles_contract(small_1v1c_problem):
    h_i, h_n, trans, eval_i, evec_i, _, _ = _exact_problem_data(small_1v1c_problem)
    spectrum, poles = rixs_krylov_scipy(
        eval_i[:1],
        evec_i[:, :1],
        h_i,
        h_n,
        trans,
        ominc=[1.0],
        eloss=[0.0, 0.2],
        gamma_c=0.2,
        gamma_f=0.1,
        nkryl=h_i.shape[0],
        workers=1,
        return_poles=True,
    )

    assert spectrum.shape == (1, 2, 1)
    assert len(poles) == 1 and len(poles[0]) == 1
    assert set(poles[0][0]) == {"npoles", "eigval", "norm", "alpha", "beta"}
    assert len(poles[0][0]["eigval"]) == 1


@pytest.mark.parallel
@pytest.mark.integration
@pytest.mark.skipif(
    "fork" not in mp.get_all_start_methods(),
    reason="stable process-pool regression uses POSIX fork",
)
def test_parallel_rixs_matches_serial(small_1v1c_problem):
    h_i, h_n, trans, eval_i, evec_i, _, _ = _exact_problem_data(small_1v1c_problem)
    common = dict(
        eval_i=eval_i[:1],
        evec_i=evec_i[:, :1],
        hmat_i=h_i,
        hmat_n=h_n,
        trans_op=trans,
        ominc=[0.8, 1.0],
        eloss=[0.0, 0.2, 0.4],
        gamma_c=0.2,
        gamma_f=0.1,
        nkryl=h_i.shape[0],
        linsys_tol=1e-11,
        linsys_maxiter=200,
        linsys_restart=max(5, h_n.shape[0]),
    )

    serial = rixs_krylov_scipy(**common, workers=1)
    parallel = rixs_krylov_scipy(
        **common, workers=2, blas_threads=1, mp_start_method="fork"
    )

    assert_allclose(parallel, serial, rtol=1e-11, atol=1e-12)
