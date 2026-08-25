"""Contract tests for the intentionally unimplemented PETSc backend stub."""

import numpy as np
import pytest

import edrixs.petsc_backend.petsc_backend as backend


def test_owns_operator_returns_false_for_non_petsc_objects():
    """Recognition must never raise when given ordinary Python/NumPy objects."""
    assert backend.owns_operator_petsc(np.eye(2)) is False
    assert backend.owns_operator_petsc(object()) is False


def test_public_symbols_are_exported():
    for name in (
        "owns_operator_petsc",
        "build_op_petsc",
        "ed_petsc",
        "xas_petsc",
        "rixs_petsc",
    ):
        assert name in backend.__all__
        assert hasattr(backend, name)


def test_build_op_stub_accepts_new_shared_builder_options(monkeypatch):
    """The shared API may pass ``use_numba`` even while PETSc is a stub."""
    monkeypatch.setattr(backend, "_petsc_module", lambda: object())
    with pytest.raises(NotImplementedError, match="build_op_petsc"):
        backend.build_op_petsc(
            None,
            None,
            object(),
            use_numba=True,
            backend_kws={"unused": True},
        )


@pytest.mark.parametrize(
    "call,operation",
    [
        (lambda: backend.ed_petsc(object()), "ed_petsc"),
        (
            lambda: backend.xas_petsc(
                np.array([0.0]), [object()], object(), [object()] * 3,
                np.array([1.0])
            ),
            "xas_petsc",
        ),
        (
            lambda: backend.rixs_petsc(
                np.array([0.0]), [object()], object(), object(),
                [object()] * 3, np.array([1.0]), np.array([0.0])
            ),
            "rixs_petsc",
        ),
    ],
)
def test_solver_stubs_remain_unimplemented(monkeypatch, call, operation):
    monkeypatch.setattr(backend, "_petsc_module", lambda: object())
    with pytest.raises(NotImplementedError, match=operation):
        call()
