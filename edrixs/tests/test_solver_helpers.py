import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse.linalg import aslinearoperator

import edrixs._solvers_helpers as helpers


def test_dense_to_sparse_u_preserves_flattening_convention():
    rng = np.random.default_rng(3)
    umat = rng.normal(size=(3, 3, 3, 3))
    umat[np.abs(umat) < 0.6] = 0

    actual = helpers._umat_dense_to_sparse(umat).toarray()

    assert_allclose(actual, umat.reshape(9, 9))


def test_sparse_siam_embedding_matches_dense_embedding():
    rng = np.random.default_rng(8)
    compact = rng.normal(size=(4, 4, 4, 4))
    dense = helpers._embed_impurity_core_umat(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )
    sparse = helpers._embed_impurity_core_umat_sparse(
        compact, v_norb=2, c_norb=2, ntot_v=6
    )

    assert_allclose(sparse.toarray(), dense.reshape(64, 64))


def test_apply_linear_combination_handles_complex_and_all_zero_coefficients():
    ops = [aslinearoperator(np.eye(2)), aslinearoperator(np.diag([2.0, 3.0]))]
    vec = np.array([1.0, -2.0])

    actual = helpers._apply_linear_combination(ops, [1j, -0.5], vec)
    expected = 1j * vec - 0.5 * (np.diag([2.0, 3.0]) @ vec)

    assert_allclose(actual, expected)
    assert_allclose(helpers._apply_linear_combination(ops, [0, 0], vec), 0)


def test_expand_broadening_accepts_scalar_or_exact_length_array():
    assert_allclose(helpers._expand_broadening(0.2, 3, "gamma"), [0.2] * 3)
    assert_allclose(helpers._expand_broadening([0.1, 0.2], 2, "gamma"), [0.1, 0.2])

    with pytest.raises(ValueError, match="shape"):
        helpers._expand_broadening([0.1], 2, "gamma")


def test_zero_rhs_short_circuits_rixs_contribution(monkeypatch):
    def forbidden_gmres(*args, **kwargs):
        raise AssertionError("GMRES must not run for a zero RHS")

    monkeypatch.setattr(helpers, "_gmres_scipy_compat", forbidden_gmres)
    record = helpers._rixs_krylov_one_contribution_scipy(
        hmat_i=aslinearoperator(np.eye(2)),
        hmat_n=aslinearoperator(np.eye(3)),
        trans_op_H=[],
        polvec_f=np.zeros(0),
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]]),
        istate=0,
        omega=1.0,
        gamma_c=0.2,
        rhs=np.zeros(3),
        skip_gs=False,
        nkryl=2,
        linsys_tol=1e-10,
        linsys_maxiter=10,
        linsys_restart=5,
    )

    assert record["norm"] == 0
    assert record["npoles"] == 1
    assert_allclose(record["alpha"], [0.0])
    assert record["beta"].size == 0


def test_gmres_failure_is_reported_with_context(monkeypatch):
    monkeypatch.setattr(
        helpers,
        "_gmres_scipy_compat",
        lambda *args, **kwargs: (np.zeros(2, dtype=complex), 7),
    )

    with pytest.raises(RuntimeError, match=r"istate=0.*omega=1.25.*info=7"):
        helpers._rixs_krylov_one_contribution_scipy(
            hmat_i=aslinearoperator(np.eye(2)),
            hmat_n=aslinearoperator(np.eye(2)),
            trans_op_H=[aslinearoperator(np.eye(2))] * 3,
            polvec_f=np.array([1.0, 0.0, 0.0]),
            eval_i=np.array([0.0]),
            evec_i=np.array([[1.0], [0.0]]),
            istate=0,
            omega=1.25,
            gamma_c=0.2,
            rhs=np.ones(2),
            skip_gs=False,
            nkryl=2,
            linsys_tol=1e-10,
            linsys_maxiter=10,
            linsys_restart=5,
        )
