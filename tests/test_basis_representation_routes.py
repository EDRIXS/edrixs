"""Tests for compact Fock-basis metadata and interchangeable realizations."""

import builtins

import numpy as np
import pytest

import edrixs._operator_builder as operator_builder
from edrixs.fock_basis import (
    FockBasis,
    FockBasisSpec,
    FockBinByN,
    build_fock_basis,
    get_fock_bin_by_N,
)
from edrixs.solvers import build_op


def _historical_integer_states(*args):
    return [
        int("".join(map(str, state)), 2)
        for state in get_fock_bin_by_N(*args)
    ]


@pytest.mark.parametrize(
    "args",
    [
        (4, 2),
        (4, 2, 2, 1),
        (2, 1, 3, 2, 2, 1),
        (3, 0, 2, 2),
    ],
)
def test_combinadic_and_explicit_match_historical_edrixs_ordering(args):
    """Both realizations must use the pre-existing EDRIXS state ordering."""
    spec = FockBasisSpec.from_args(*args)
    explicit = build_fock_basis(spec, method="explicit")
    combinadic = build_fock_basis(spec, method="combinadic")
    expected = _historical_integer_states(*args)

    assert isinstance(explicit, FockBasis)
    assert isinstance(combinadic, FockBinByN)
    assert explicit.basis_int == expected
    assert len(explicit) == len(combinadic) == len(expected)

    for index, state in enumerate(expected):
        assert explicit.decode(index) == state
        assert combinadic.decode(index) == state
        assert explicit.encode(state) == index
        assert combinadic.encode(state) == index


def test_basis_spec_is_compact_sector_metadata():
    """A basis specification stores shell structure without materializing states."""
    spec = FockBasisSpec.from_args(4, 2, 2, 1)

    assert spec.shapes == ((4, 2), (2, 1))
    assert spec.norbs == 6
    assert len(spec) == 12
    assert not hasattr(spec, "basis_int")
    assert not hasattr(spec, "lookup")


def test_build_fock_basis_defaults_to_combinadic_and_can_select_explicit():
    spec = FockBasisSpec.from_args(4, 2)

    assert isinstance(build_fock_basis(spec), FockBinByN)
    assert isinstance(build_fock_basis(spec, method="explicit"), FockBasis)


def test_realized_arbitrary_explicit_basis_is_preserved():
    """SciPy-compatible arbitrary explicit subspaces remain valid inputs."""
    arbitrary = FockBasis([0b1100, 0b1001, 0b0011], norbs=4)

    assert build_fock_basis(arbitrary, method="explicit") is arbitrary
    assert build_fock_basis(arbitrary, method="combinadic") is arbitrary


def test_invalid_basis_method_is_reported():
    with pytest.raises(ValueError, match="combinadic.*explicit"):
        build_fock_basis(FockBasisSpec.from_args(4, 2), method="unknown")


def test_numba_is_only_required_when_explicitly_requested(monkeypatch):
    """The default operator route must work when importing numba is impossible."""
    spec = FockBasisSpec.from_args(2, 1)
    emat = np.array([[1.0, 0.25], [0.25, -0.5]], dtype=complex)

    monkeypatch.setattr(operator_builder, "_NUMBA_KERNELS", None)
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "numba" or name.startswith("numba."):
            raise ImportError("numba intentionally unavailable in this test")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    operator = build_op(
        emat,
        None,
        spec,
        backend="scipy",
        basis_method="combinadic",
        use_numba=False,
    )
    np.testing.assert_allclose(operator.toarray(), emat)

    with pytest.raises(ImportError, match="use_numba=True requires numba"):
        build_op(
            emat,
            None,
            spec,
            backend="scipy",
            basis_method="combinadic",
            use_numba=True,
        )
