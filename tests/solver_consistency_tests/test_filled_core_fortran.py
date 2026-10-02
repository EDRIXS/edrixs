"""Modern filled-core and historical legacy Fortran workflows."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

import edrixs

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(importlib.util.find_spec('edrixs.fedrixs') is None,
                       reason='requires compiled Fortran solvers'),
    pytest.mark.filterwarnings('ignore:.*is deprecated.*:DeprecationWarning'),
]


def _case(kind):
    # Break spin degeneracy so one-state spectra do not depend on eigenvectors
    # chosen inside a degenerate ground-state subspace.
    common = dict(slater=([0.4] * 30, [0.5] * 30), c_soc=0.11)
    if kind == '1v1c':
        return dict(shell_name=('s', 'p'), v_noccu=1, shell_level=(0.2, -3),
                    v_othermat=np.array([[-0.2, 0.03j], [-0.03j, 0.2]]), **common)
    if kind == '2v1c':
        return dict(shell_name=('s', 's', 'p'), v_tot_noccu=1,
                    shell_level=(0.2, 0.7, -3),
                    v1_othermat=np.diag([-0.2, 0.2]),
                    hopping_v1v2=np.array([[0.1, 0.03j], [0.02j, 0.08]]), **common)
    return dict(shell_name=('s', 'p'), nbath=1, v_noccu=1, c_level=-3,
                imp_mat=np.diag([-0.2, 0.2]), bath_level=np.array([[0.6, 0.8]]),
                hyb=np.array([[0.1+0.02j, 0.08-0.03j]]), **common)


@pytest.mark.parametrize('kind', ['1v1c', '2v1c', 'siam'])
@pytest.mark.parametrize('route', ['modern', 'legacy'])
def test_fortran_ed_xas_rixs_matches_its_core_convention(kind, route, tmp_path, monkeypatch):
    from mpi4py import MPI
    monkeypatch.chdir(tmp_path)
    kwargs = _case(kind)
    problem = list(getattr(edrixs, 'model_' + kind)(**kwargs))
    nv, nocc = problem[2].shapes[0]
    if route == 'legacy':
        # Historical native initial states omit core interactions. Retaining
        # the filled core's one-body energy gives the same legacy reference.
        valence_u = problem[1][:nv, :nv, :nv, :nv].copy()
        problem[1] = np.zeros_like(problem[1])
        problem[1][:nv, :nv, :nv, :nv] = valence_u
    dense_i, dense_n, dense_t = edrixs.get_ops(*problem[:7], backend='dense')
    exact_e, exact_v = np.linalg.eigh(dense_i)
    center = np.median(np.linalg.eigvalsh(dense_n)) - exact_e[0]
    ominc = np.array([center - 0.15, center + 0.2])
    loss = np.linspace(-0.1, 1.0, 5)
    xargs = dict(gamma_c=0.25, thin=0.7, phi=0.2, pol_type=[('isotropic', 0)])
    rargs = dict(gamma_c=0.25, gamma_f=0.1, thin=0.7, thout=1.1, phi=0.2,
                 pol_type=[('linear', 0.1, 'linear', -0.2)])
    exact_x = edrixs.xas(exact_e[:1], exact_v[:, :1], dense_n, dense_t, ominc,
                         backend='dense', **xargs)
    exact_r = edrixs.rixs(exact_e[:1], exact_v[:, :1], dense_i, dense_n, dense_t,
                          ominc, loss, backend='dense', **rargs)
    if route == 'modern':
        hi, hn, trans = edrixs.get_ops(*problem[:7], backend='fortran')
        energies, vectors = edrixs.ed(hi, num_evals=1,
                                      backend_kws={'ed_solver': 0, 'nvector': 1})
        absorption = edrixs.xas(energies, vectors, hn, trans, ominc,
                                backend_kws={'nkryl': len(dense_n)}, **xargs)
        scattering = edrixs.rixs(
            energies, vectors, hi, hn, trans, ominc, loss,
            backend_kws={'nkryl': len(dense_n), 'linsys_tol': 1e-12}, **rargs)
    else:
        legacy_ed = getattr(edrixs, 'ed_' + kind + '_fort')
        result = legacy_ed(MPI.COMM_WORLD, **kwargs, do_ed=2 if kind == 'siam' else False)
        assert all(item is None for item in result)
        # Native legacy files retain the full uncontracted Coulomb input and
        # redistribute only the filled core's one-body energy over valence.
        original = getattr(edrixs, 'model_' + kind)(**kwargs)
        expected_h = np.zeros_like(original[0])
        expected_h[:nv, :nv] = original[0][:nv, :nv]
        expected_h[:nv, :nv] += np.eye(nv) * np.trace(original[0][nv:, nv:]) / nocc
        edrixs.write_emat(expected_h, 'expected_h.in')
        edrixs.write_umat(original[1], 'expected_u.in')
        for actual, expected in [('hopping_i.in', 'expected_h.in'),
                                 ('coulomb_i.in', 'expected_u.in')]:
            assert_allclose(np.loadtxt(actual, skiprows=1),
                            np.loadtxt(expected, skiprows=1), atol=1e-13)
        staged = {name: Path(name).read_bytes() for name in (
            'hopping_i.in', 'hopping_n.in', 'coulomb_i.in', 'coulomb_n.in', 'fock_i.in')}
        legacy = legacy_ed(MPI.COMM_WORLD, **kwargs, ed_solver=0, neval=1,
                           nvector=1, idump=True, do_ed=1)
        for name, contents in staged.items():
            assert Path(name).read_bytes() == contents
        energies, denmat = legacy[:2]
        assert denmat.shape == (1, nv, nv)
        assert_allclose(np.trace(denmat[0]), nocc, atol=1e-12)
        spectrum_args = dict(shell_name=kwargs['shell_name'])
        if kind == 'siam':
            spectrum_args['nbath'] = 1
        absorption, _ = getattr(edrixs, 'xas_' + kind + '_fort')(
            MPI.COMM_WORLD, ominc=ominc, nkryl=len(dense_n), **spectrum_args, **xargs)
        scattering, _ = getattr(edrixs, 'rixs_' + kind + '_fort')(
            MPI.COMM_WORLD, ominc=ominc, eloss=loss, nkryl=len(dense_n),
            linsys_tol=1e-12, **spectrum_args, **rargs)
        assert Path('fock_f.in').read_bytes() == staged['fock_i.in']
    assert_allclose(energies, exact_e[:1], atol=2e-9)
    assert_allclose(absorption, exact_x, rtol=2e-7, atol=2e-9)
    assert_allclose(scattering, exact_r, rtol=2e-6, atol=2e-8)


def test_legacy_siam_occupancy_search_keeps_historical_core_convention(tmp_path, monkeypatch):
    from mpi4py import MPI
    monkeypatch.chdir(tmp_path)
    # Legacy searches omit the core potential and constant. These parameters
    # fill both impurity and bath (N=4); the modern filled-core model has N=2.
    kwargs = dict(shell_name=('s', 's'), nbath=1, c_level=-3,
                  imp_mat=np.diag([-0.8, -0.7]), bath_level=np.array([[-0.5, -0.4]]),
                  slater=([0, 0.8, 0, 0.6], [0, 0.8, 0, 0.6]))
    exact = {}
    for nocc in range(5):
        problem = edrixs.model_siam(**kwargs, v_noccu=nocc)
        hi = edrixs.build_op(
            problem[0][:4, :4], problem[1][:4, :4, :4, :4],
            edrixs.FockBasisSpec.from_args(4, nocc), backend='dense')
        exact[nocc] = np.linalg.eigvalsh(hi)[0]
    assert min(exact, key=exact.get) == 4
    energies, density, nocc = edrixs.ed_siam_fort(
        MPI.COMM_WORLD, **kwargs, do_ed=0, ed_solver=0, neval=1, nvector=1)
    assert nocc == 4
    assert_allclose(energies, [exact[nocc] + 2 * kwargs['c_level']], atol=2e-9)
    assert_allclose(np.trace(density[0]), nocc, atol=1e-12)
    for filename in ('search_gs.log', 'search_result.dat'):
        for sector, energy, _ in np.loadtxt(filename, ndmin=2):
            assert_allclose(energy, exact[int(sector)], atol=2e-9)
    staged = Path('hopping_i.in').read_bytes()
    edrixs.ed_siam_fort(MPI.COMM_WORLD, **kwargs, v_noccu=nocc, do_ed=2)
    assert Path('hopping_i.in').read_bytes() == staged


@pytest.mark.parametrize('kind', ['1v1c', '2v1c', 'siam'])
def test_legacy_file_only_modes_do_not_write_on_nonroot(kind, tmp_path, monkeypatch):
    import edrixs.solvers as solvers
    import edrixs._solvers_helpers as helpers
    from edrixs.fortran_backend import fortran_backend

    class Nonroot:
        def Get_rank(self):
            return 1

        def Get_size(self):
            return 2

        def py2f(self):
            return 0

    def forbidden(*args, **kwargs):
        raise AssertionError('non-root wrote native input files')

    monkeypatch.chdir(tmp_path)
    for module in (solvers, helpers, fortran_backend):
        for name in ('write_emat', 'write_umat', 'write_config'):
            monkeypatch.setattr(module, name, forbidden)
    monkeypatch.setattr(solvers, 'write_fock_dec_by_N', forbidden)
    monkeypatch.setattr(helpers, 'write_fock_dec_by_N', forbidden)
    result = getattr(edrixs, 'ed_' + kind + '_fort')(
        Nonroot(), **_case(kind), do_ed=2 if kind == 'siam' else False)
    assert all(item is None for item in result)
    assert not list(tmp_path.iterdir())
