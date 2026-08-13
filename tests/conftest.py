"""Small shared physical model used by unit and consistency tests."""

import numpy as np
import pytest

from edrixs.models import model_1v1c
from edrixs.solvers import get_ops


@pytest.fixture(scope="module")
def small_1v1c_kwargs():
    """Return a tiny dipole-allowed model that keeps workflows fast."""
    return {
        "shell_name": ("p", "s"),
        "shell_level": (0.2, -4.0),
        "v_soc": (0.05, 0.08),
        "v_noccu": 1,
        "slater": ([0.7], [0.9]),
        "v_cfmat": np.diag(np.linspace(-0.12, 0.12, 6)),
    }


@pytest.fixture(scope="module")
def small_1v1c_problem(small_1v1c_kwargs):
    """Build the shared orbital-space model before selecting a backend."""
    return model_1v1c(**small_1v1c_kwargs, sparse_U=False)


@pytest.fixture(scope="module")
def small_1v1c_operators(small_1v1c_problem):
    """Build SciPy Hamiltonians and transition operators for public-API tests."""
    return get_ops(*small_1v1c_problem, backend="scipy")
