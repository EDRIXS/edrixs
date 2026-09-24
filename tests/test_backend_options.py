"""Contract tests for staged ``backend_kws`` validation."""

import inspect

import numpy as np
import pytest

from edrixs.dense_backend import dense_backend
from edrixs.fortran_backend import fortran_backend
from edrixs.fortran_backend.options import OPTIONS as FORTRAN_OPTIONS
from edrixs.dense_backend.options import OPTIONS as DENSE_OPTIONS
from edrixs.petsc_backend.options import OPTIONS as PETSC_OPTIONS
from edrixs.scipy_backend import scipy_backend
from edrixs.scipy_backend.options import OPTIONS as SCIPY_OPTIONS
from edrixs.solvers import rixs_1v1c_fort, rixs_2v1c_fort, rixs_siam_fort


def test_option_registries_cover_each_staged_operation():
    assert set(DENSE_OPTIONS) == {'build_op', 'get_ops', 'ed', 'xas', 'rixs'}
    assert set(SCIPY_OPTIONS) == {'build_op', 'get_ops', 'ed', 'xas', 'rixs'}
    assert set(FORTRAN_OPTIONS) == {'get_ops', 'ed', 'xas', 'rixs'}
    assert set(PETSC_OPTIONS) == {'build_op', 'get_ops', 'ed', 'xas', 'rixs'}
    assert all(not options for options in PETSC_OPTIONS.values())
    assert not DENSE_OPTIONS['ed']
    assert not DENSE_OPTIONS['xas']
    assert not DENSE_OPTIONS['rixs']


def test_scipy_unknown_option_suggests_close_name():
    with pytest.raises(TypeError, match=r"maxter.*maxiter"):
        scipy_backend.ed_scipy(
            np.eye(2), backend_kws={'maxter': 10}
        )


def test_option_from_wrong_operation_reports_valid_operation():
    with pytest.raises(TypeError, match=r"blocksize.*scipy\.ed"):
        scipy_backend.xas_scipy(
            None, None, None, None, None,
            backend_kws={'blocksize': 2},
        )


def test_backend_options_require_string_keys():
    with pytest.raises(TypeError, match='keys must be strings'):
        scipy_backend.ed_scipy(np.eye(2), backend_kws={1: 2})


@pytest.mark.parametrize(
    'backend_kws, error',
    [
        ({'maxiter': 0}, ValueError),
        ({'tol': 'small'}, TypeError),
        ({'tol': float('nan')}, ValueError),
        ({'suppress_lobpcg_warnings': 1}, TypeError),
    ],
)
def test_scipy_ed_option_values_are_validated(backend_kws, error):
    with pytest.raises(error):
        scipy_backend.ed_scipy(np.eye(2), backend_kws=backend_kws)


def test_dense_solver_rejects_iterative_scipy_options():
    """Dense eigensums do not expose Krylov or GMRES controls."""
    with pytest.raises(TypeError, match=r"dense\.ed.*tol"):
        dense_backend.ed_dense(
            np.eye(2), backend_kws={'tol': 1e-10}
        )


def test_fortran_staged_rixs_uses_standardized_iteration_name():
    assert 'linsys_maxiter' in FORTRAN_OPTIONS['rixs']
    assert 'linsys_max' not in FORTRAN_OPTIONS['rixs']

    with pytest.raises(TypeError, match=r"linsys_max.*linsys_maxiter"):
        fortran_backend.rixs_fortran(
            [], None, None, None, None, [], [],
            backend_kws={'linsys_max': 10},
        )


def test_legacy_fortran_rixs_signatures_are_unchanged():
    for function in (rixs_1v1c_fort, rixs_2v1c_fort, rixs_siam_fort):
        parameters = inspect.signature(function).parameters
        assert 'linsys_max' in parameters
        assert 'linsys_maxiter' not in parameters
