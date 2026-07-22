import numpy as np
import pytest

from edrixs.solvers import setup_1v1c


@pytest.fixture(scope="module")
def small_1v1c_kwargs():
    # A small dipole-allowed p <- s problem with nonzero one- and two-body terms.
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
    return setup_1v1c(**small_1v1c_kwargs, sparse_U=False)
