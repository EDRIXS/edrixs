import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator

import edrixs._solvers_helpers as helpers
from edrixs.solvers import rixs_krylov_scipy, xas_krylov_scipy


def _valid_xas_kwargs():
    return dict(
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]], dtype=complex),
        hmat_n=np.eye(3),
        trans_op=[np.zeros((3, 2), dtype=complex) for _ in range(3)],
        ominc=np.array([0.8, 1.0]),
        gamma_c=0.2,
        nkryl=3,
    )


@pytest.mark.parametrize(
    "update, message",
    [
        ({"eval_i": np.array([[0.0]])}, "one-dimensional"),
        ({"evec_i": np.array([1.0, 0.0])}, "two-dimensional"),
        ({"evec_i": np.eye(2, dtype=complex)}, "columns"),
        ({"hmat_n": np.zeros((3, 2))}, "square"),
        ({"trans_op": [np.zeros((3, 2)) for _ in range(2)]}, "len\\(trans_op\\)"),
        (
            {"trans_op": [np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((2, 2))]},
            "expected",
        ),
        ({"scatter_axis": np.eye(2)}, "scatter_axis"),
        ({"nkryl": 0}, "positive"),
        ({"gamma_c": [0.2]}, "shape"),
    ],
)
def test_xas_public_validation(update, message):
    kwargs = _valid_xas_kwargs()
    kwargs.update(update)
    with pytest.raises(ValueError, match=message):
        xas_krylov_scipy(**kwargs)


def test_xas_rejects_unknown_polarization():
    kwargs = _valid_xas_kwargs()
    kwargs["pol_type"] = [("diagonal", 0.0)]
    with pytest.raises(ValueError, match="Unknown XAS polarization"):
        xas_krylov_scipy(**kwargs)


def _valid_rixs_kwargs():
    transitions = [np.zeros((3, 2), dtype=complex) for _ in range(3)]
    transitions[0][0, 0] = 1.0
    return dict(
        eval_i=np.array([0.0]),
        evec_i=np.array([[1.0], [0.0]], dtype=complex),
        hmat_i=np.eye(2),
        hmat_n=np.eye(3),
        trans_op=transitions,
        ominc=np.array([1.0]),
        eloss=np.array([0.0, 0.2]),
        gamma_c=0.2,
        gamma_f=0.1,
        nkryl=2,
        linsys_maxiter=20,
        linsys_restart=5,
        workers=1,
    )


@pytest.mark.parametrize(
    "update, message",
    [
        ({"hmat_i": np.zeros((2, 3))}, "hmat_i must be square"),
        ({"hmat_n": np.zeros((3, 2))}, "hmat_n must be square"),
        ({"eval_i": np.array([[0.0]])}, "one-dimensional"),
        ({"evec_i": np.array([1.0, 0.0])}, "two-dimensional"),
        (
            {"eval_i": np.array([0.0, 0.2])},
            "len\\(eval_i\\) must equal",
        ),
        ({"evec_i": np.ones((3, 1))}, "evec_i.shape\\[0\\]"),
        ({"trans_op": [np.zeros((3, 2)) for _ in range(2)]}, "len\\(trans_op\\)"),
        (
            {"trans_op": [np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((2, 2))]},
            "expected",
        ),
        ({"scatter_axis": np.eye(2)}, "scatter_axis"),
        ({"nkryl": 0}, "nkryl"),
        ({"linsys_maxiter": 0}, "linsys_maxiter"),
        ({"linsys_restart": 0}, "linsys_restart"),
        ({"blas_threads": 0}, "blas_threads"),
        ({"workers": 0}, "workers"),
        ({"gamma_c": [0.1, 0.2]}, "gamma_c.*shape"),
        ({"gamma_f": [0.1]}, "gamma_f.*shape"),
    ],
)
def test_rixs_public_validation(update, message):
    kwargs = _valid_rixs_kwargs()
    kwargs.update(update)
    with pytest.raises(ValueError, match=message):
        rixs_krylov_scipy(**kwargs)


def test_rixs_requires_transition_adjoint_action():
    # A rectangular LinearOperator with matvec only cannot supply T.H @ x.
    no_adjoint = LinearOperator(
        shape=(3, 2),
        matvec=lambda x: np.zeros(3, dtype=complex),
        dtype=np.complex128,
    )
    kwargs = _valid_rixs_kwargs()
    kwargs["trans_op"] = [no_adjoint, no_adjoint, no_adjoint]

    with pytest.raises(TypeError, match="adjoint action"):
        rixs_krylov_scipy(**kwargs)


def test_process_pool_rejects_nonserializable_operator():
    nonserializable = LinearOperator(
        shape=(2, 2),
        matvec=lambda x: x,
        dtype=np.complex128,
    )

    with pytest.raises(TypeError, match="not serializable"):
        helpers._operator_for_process_pool(nonserializable, "custom_op")


def test_rixs_worker_job_requires_initializer(monkeypatch):
    monkeypatch.setattr(helpers, "_RIXS_WORKER_STATE", None)
    with pytest.raises(RuntimeError, match="not initialized"):
        helpers._rixs_pool_job((0, 0, 0))
