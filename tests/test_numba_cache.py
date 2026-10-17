"""Verify that operator kernels reuse compiled code across Python processes."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.skipif(importlib.util.find_spec("numba") is None,
                    reason="Numba is required for disk caching")
def test_operator_kernels_load_from_disk_in_fresh_process(tmp_path):
    script = """
import sys
import numpy as np
from edrixs import build_op
from edrixs.fock_basis import FockBasisSpec
from edrixs._operator_builder import _get_numba_kernels

for method in ('combinadic', 'explicit'):
    operator = build_op(np.diag([1., 2.]), None,
                        FockBasisSpec.from_args(2, 1),
                        backend='scipy', basis_method=method)()
    np.testing.assert_allclose(np.sort(operator.diagonal()), [1., 2.])
for kernel in _get_numba_kernels():
    assert sum(kernel.stats.cache_hits.values()) == int(sys.argv[1])
    assert sum(kernel.stats.cache_misses.values()) == 1 - int(sys.argv[1])
"""
    env = dict(os.environ, NUMBA_CACHE_DIR=str(tmp_path / "cache"))
    root = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [root, env.get("PYTHONPATH")]))
    for expected_hits in (0, 1):
        result = subprocess.run(
            [sys.executable, "-c", script, str(expected_hits)],
            env=env, cwd=root, capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0, result.stdout + result.stderr
