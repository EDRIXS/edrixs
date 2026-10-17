"""The public solver interface consumes labeled, deferred matrix factories."""

from functools import partial

import numpy as np
import pytest

from edrixs import FockBasisSpec, build_op, ed, get_ops, rixs, xas
from edrixs import solvers
from edrixs.fortran_backend.fortran_backend import FortranDiskOperator


def labeled(builder, backend='scipy'):
    factory = partial(builder)
    factory.backend = backend
    return factory


@pytest.mark.parametrize('backend', ['dense', 'scipy', 'petsc'])
def test_build_op_defers_backend_construction(backend, monkeypatch):
    calls = []
    matrix = object()
    module = getattr(solvers, backend + '_backend')

    def builder(*args, **kwargs):
        calls.append((args, kwargs))
        return matrix

    monkeypatch.setattr(module, 'build_op_' + backend, builder)
    factory = build_op(np.eye(2), None, FockBasisSpec.from_args(2, 1),
                       backend=backend)
    assert isinstance(factory, partial)
    assert factory.backend == backend
    assert calls == []
    assert factory() is matrix
    assert len(calls) == 1


@pytest.mark.parametrize('backend', ['dense', 'scipy'])
def test_factory_builds_fresh_matrices_and_retains_coefficient_references(backend):
    coefficients = np.diag([1., 2.])
    factory = build_op(coefficients, None, FockBasisSpec.from_args(2, 1),
                       backend=backend, use_numba=False)
    first = factory()
    coefficients[0, 0] = 3.
    second = factory()
    assert first is not second
    assert sorted(first.diagonal()) == [1., 2.]
    assert sorted(second.diagonal()) == [2., 3.]


@pytest.mark.parametrize('solver,bad_kind', [
    (solver, bad_kind)
    for solver in ['ed', 'xas', 'rixs']
    for bad_kind in ['missing', 'unknown', 'mixed']
    if solver != 'ed' or bad_kind != 'mixed'
])
def test_labels_are_validated_before_any_construction(solver, bad_kind):
    calls = []
    good = labeled(lambda: calls.append(1))
    bad = (object() if bad_kind == 'missing' else
           labeled(lambda: calls.append(2), 'BAD' if bad_kind == 'unknown' else 'dense'))
    expected = TypeError if bad_kind == 'missing' else ValueError
    with pytest.raises(expected):
        if solver == 'ed':
            ed(bad)
        elif solver == 'xas':
            xas(None, None, good, [good, bad, good], None)
        else:
            rixs(None, None, good, good, [good, bad, good], None, None)
    assert calls == []


@pytest.mark.parametrize('solver', ['ed', 'xas', 'rixs'])
def test_materialize_each_distinct_factory_once_per_solver_call(solver, monkeypatch):
    calls = []
    matrix = object()
    factory = labeled(lambda: calls.append(1) or matrix)
    result = object()

    def backend(*args, **kwargs):
        if solver == 'ed':
            assert args[0] is matrix
        elif solver == 'xas':
            assert args[2] is matrix
            assert all(op is matrix for op in args[3])
        else:
            assert args[2] is matrix and args[3] is matrix
            assert all(op is matrix for op in args[4])
        return result

    monkeypatch.setattr(solvers.scipy_backend, solver + '_scipy', backend)

    def run():
        if solver == 'ed':
            return ed(factory)
        if solver == 'xas':
            return xas(None, None, factory, [factory] * 3, None)
        return rixs(None, None, factory, factory, [factory] * 3, None, None)
    assert run() is result
    assert calls == [1]
    assert run() is result
    assert calls == [1, 1]


def test_distinct_factories_build_independently(monkeypatch):
    matrices = [object(), object()]
    calls = []
    factories = [labeled(lambda i=i: calls.append(i) or matrices[i]) for i in range(2)]

    def backend(energies, states, hamiltonian, transitions, grid, **kwargs):
        assert hamiltonian is matrices[0]
        assert transitions == [matrices[1], matrices[0], matrices[1]]
    monkeypatch.setattr(solvers.scipy_backend, 'xas_scipy', backend)
    xas(None, None, factories[0], [factories[1], factories[0], factories[1]], None)
    assert calls == [0, 1]


@pytest.mark.parametrize('solver', ['ed', 'xas', 'rixs'])
def test_fortran_handles_pass_through_without_calling(solver, monkeypatch):
    handle = FortranDiskOperator()
    assert handle.backend == 'fortran'
    result = object()

    def backend(*args, **kwargs):
        index = 0 if solver == 'ed' else 2
        assert args[index] is handle
        return result
    monkeypatch.setattr(solvers.fortran_backend, solver + '_fortran', backend)
    if solver == 'ed':
        actual = ed(handle)
    elif solver == 'xas':
        actual = xas(None, None, handle, [handle] * 3, None)
    else:
        actual = rixs(None, None, handle, handle, [handle] * 3, None, None)
    assert actual is result


def test_solver_backend_keyword_removed():
    factory = labeled(lambda: np.eye(2))
    for solver, args in [(ed, [factory]),
                         (xas, [None, None, factory, [factory] * 3, None]),
                         (rixs, [None, None, factory, factory, [factory] * 3, None, None])]:
        with pytest.raises(TypeError, match="unexpected keyword argument 'backend'"):
            solver(*args, backend='scipy')


def test_backend_options_validation_occurs_when_building():
    factory = build_op(np.eye(2), None, FockBasisSpec.from_args(2, 1),
                       backend_kws={'unknown_option': 1})
    with pytest.raises(TypeError, match='Unknown backend_kws'):
        factory()


def test_ed_rejects_nonmapping_solver_options():
    factory = labeled(lambda: np.eye(2))
    with pytest.raises(TypeError, match='mapping'):
        ed(factory, backend_kws=[('tol', 1e-8)])


def test_get_ops_returns_labeled_factories(small_1v1c_problem):
    hi, hn, transitions = get_ops(*small_1v1c_problem[:7], use_numba=False)
    assert all(isinstance(op, partial) and op.backend == 'scipy'
               for op in [hi, hn, *transitions])
    values, vectors = ed(hi, num_evals=1)
    assert len(values) == 1
    assert vectors.shape[1] == 1
