"""Exercise example workflows against the public solver API on a small model."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from edrixs.models import model_1v1c
from edrixs.solvers import get_ops


EXAMPLES = Path(__file__).resolve().parents[1] / 'examples' / 'scipy'


@pytest.mark.parametrize('script', sorted(EXAMPLES.glob('*.py')), ids=lambda p: p.stem)
def test_scipy_example_workflow(script, tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(script.stem, script)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)

    # Keep enough states for every example's ED block, with a nonzero shift
    # to check residuals on the physical energy reference.
    problem = model_1v1c(
        ('d', 'p'), v_noccu=3, shell_level=(0.0, -2.0),
        v_soc=(0.1, 0.1), c_soc=0.2,
        v_cfmat=np.diag(np.linspace(-0.3, 0.3, 10)), sparse_U=True,
    )
    for name in ('model_1v1c', 'model_2v1c', 'model_siam'):
        if hasattr(example, name):
            monkeypatch.setattr(example, name, lambda *args, **kwargs: problem)
    monkeypatch.setattr(example, 'get_ops',
                        lambda *args, **kwargs: get_ops(*args, **kwargs, use_numba=False))

    # Shrink only the spectral grids, retaining the real solvers and file output.
    numpy = SimpleNamespace(**{name: getattr(np, name) for name in dir(np)})
    numpy.linspace = lambda start, stop, num: np.linspace(start, stop, min(num, 12))
    monkeypatch.setattr(example, 'np', numpy)
    example.run(tmp_path)

    energies = np.loadtxt(tmp_path / 'eval_i.dat')
    assert np.isfinite(energies).all()
    assert np.max(energies[:, 2]) < 1e-8
    with np.load(next(tmp_path.glob('*_results.npz'))) as results:
        for name in ('xas', 'rixs'):
            if name == 'rixs' and 'nio_aim' in script.stem:
                continue
            spectrum = results[name]
            assert spectrum.size
            assert np.isfinite(spectrum).all()
            assert np.max(spectrum) > 0
    assert list(tmp_path.glob('*.pdf'))
