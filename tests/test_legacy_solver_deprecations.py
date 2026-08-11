"""Deprecation behavior and compatibility for legacy solver workflows."""

import inspect

import pytest

from edrixs.solvers import ed_1v1c_fort, ed_1v1c_py, rixs_1v1c_py


def test_legacy_python_solver_warns_before_validation():
    """Legacy dense-Python workflows emit a library deprecation warning."""
    with pytest.warns(DeprecationWarning, match="ed_1v1c_py is deprecated"):
        with pytest.raises(Exception, match="NOT supported"):
            ed_1v1c_py(("bad", "s"))


def test_legacy_fortran_solver_warns_before_validation():
    """Legacy Fortran workflows emit the same deprecation category."""
    with pytest.warns(DeprecationWarning, match="ed_1v1c_fort is deprecated"):
        with pytest.raises(Exception, match="NOT supported"):
            ed_1v1c_fort(object(), ("bad", "s"))


def test_legacy_rixs_retains_skip_gs_signature():
    """The deprecated dense RIXS API retains its historical elastic flag."""
    parameters = inspect.signature(rixs_1v1c_py).parameters
    assert parameters["skip_gs"].default is False
