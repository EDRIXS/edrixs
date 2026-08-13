"""Backend-neutral solver interface and deprecated legacy workflows.

Physical model construction lives in :mod:`edrixs.models`. Numerical operator
construction and spectral solvers are implemented in ``*_backend.py`` modules.
"""

from __future__ import annotations

import warnings

import numpy as np
import scipy

from .iostream import (
    write_tensor, write_emat, write_umat, write_config, read_poles_from_file
)
from .angular_momentum import (
    get_sx, get_sy, get_sz, get_lx, get_ly, get_lz,
    rmat_to_euler, get_wigner_dmat,
)
from .photon_transition import (
    get_trans_oper, quadrupole_polvec, dipole_polvec_xas,
    dipole_polvec_rixs, unit_wavevector,
)
from .coulomb_utensor import get_umat_slater
from .manybody_operator import two_fermion, four_fermion
from .fock_basis import get_fock_bin_by_N, write_fock_dec_by_N
from .basis_transform import cb_op2, tmat_r2c, cb_op
from .utils import info_atomic_shell, slater_integrals_name, boltz_dist
from .rixs_utils import scattering_mat
from .plot_spectrum import get_spectra_from_poles, merge_pole_dicts
from .soc import atom_hsoc
from .petsc_backend import petsc_backend
from .scipy_backend import scipy_backend
from ._solvers_helpers import (
    _ed_1or2_valence_1core,
    _infer_backend,
    _rixs_1or2_valence_1core,
    _xas_1or2_valence_1core,
)

__all__ = [
    # Backend-neutral staged interface.
    'build_op', 'get_ops', 'ed', 'xas', 'rixs',

    # Legacy dense-Python interface.
    'ed_1v1c_py', 'xas_1v1c_py', 'rixs_1v1c_py',

    # Legacy Fortran interface.
    'ed_1v1c_fort', 'xas_1v1c_fort', 'rixs_1v1c_fort',
    'ed_2v1c_fort', 'xas_2v1c_fort', 'rixs_2v1c_fort',
    'ed_siam_fort', 'xas_siam_fort', 'rixs_siam_fort',
]


# -----------------------------------------------------------------------------
# Backend-neutral staged interface
# -----------------------------------------------------------------------------


def build_op(emat, umat, lb, rb=None, *, backend='scipy', backend_kws=None):
    """
    Build a many-body operator with the selected backend.

    Parameters
    ----------
    emat : array-like or None
        Coefficients of the one-body part. Pass ``None`` when the operator has
        no one-body contribution.
    umat : array-like, sparse matrix, or None
        Coefficients of the two-body part. Pass ``None`` when the operator has
        no two-body contribution.
    lb : FockBasis
        Basis for the output (left) many-body space.
    rb : FockBasis or None, optional
        Basis for the input (right) many-body space. When omitted, ``lb`` is
        used for both sides.
    backend : str, optional
        Backend name. The default is ``'scipy'``.
    backend_kws : mapping, optional
        Backend-specific construction options. For the SciPy backend this
        includes ``tol``.

    Returns
    -------
    operator
        Many-body operator in the representation owned by ``backend``.
    """
    match backend:
        case 'scipy':
            build_op_backend = scipy_backend.build_op_scipy
        case 'dense':
            build_op_backend = scipy_backend.build_op_dense
        case 'petsc':
            build_op_backend = petsc_backend.build_op_petsc
        case _:
            raise ValueError(
                "Unknown backend {!r}; expected 'scipy', 'dense', or "
                "'petsc'".format(backend)
            )

    return build_op_backend(
        emat,
        umat,
        lb,
        rb,
        backend_kws=backend_kws,
    )


def get_ops(
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, *,
    backend='scipy', backend_kws=None,
):
    """
    Build initial/intermediate Hamiltonians and transition operators.

    Parameters
    ----------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat
        Backend-neutral problem definition returned by :mod:`edrixs.models`
        model functions.
    backend : {'scipy', 'dense'}, optional
        Backend used for the returned operators. The default is ``'scipy'``.
    backend_kws : mapping, optional
        Backend-specific operator-construction options. For the SciPy and
        dense compatibility backends this includes ``tol``.

    Returns
    -------
    hmat_i, hmat_n, trans_ops
        Initial/final Hamiltonian, intermediate Hamiltonian, and transition
        operators for the selected backend.
    """
    hmat_i = build_op(
        emat_i, umat_i, basis_i, backend=backend,
        backend_kws=backend_kws,
    )
    hmat_n = build_op(
        emat_n, umat_n, basis_n, backend=backend,
        backend_kws=backend_kws,
    )

    trans_mat = np.asarray(trans_mat)
    if trans_mat.ndim != 3:
        raise ValueError("trans_mat must be a three-dimensional array")
    trans_ops = [
        build_op(
            component,
            None,
            basis_n,
            basis_i,
            backend=backend,
            backend_kws=backend_kws,
        )
        for component in trans_mat
    ]
    return hmat_i, hmat_n, trans_ops


def ed(hmat_i, num_evals=1, *, backend=None, backend_kws=None):
    """
    Compute low-energy initial states through a numerical backend.

    Parameters
    ----------
    hmat_i : backend-owned operator
        Initial/final Hamiltonian.
    num_evals : int, optional
        Number of lowest eigenpairs to return.
    backend : str or None, optional
        Backend name. When omitted, infer it from ``hmat_i``.
    backend_kws : mapping, optional
        Backend-specific eigensolver options.

    Returns
    -------
    eigenvalues, eigenvectors
        Lowest retained eigenpairs.
    """
    backend_name = backend if backend is not None else _infer_backend(hmat_i)
    match backend_name:
        case 'scipy':
            ed_backend = scipy_backend.ed_scipy
        case 'dense':
            ed_backend = scipy_backend.ed_dense
        case 'petsc':
            ed_backend = petsc_backend.ed_petsc
        case _:
            raise ValueError(
                "Unknown backend {!r}; expected 'scipy', 'dense', or "
                "'petsc'".format(backend_name)
            )

    return ed_backend(
        hmat_i,
        num_evals=num_evals,
        backend_kws=backend_kws,
    )


def xas(eval_i, evec_i, hmat_n, trans_op, ominc, *,
        gamma_c=0.1, thin=1.0, phi=0.0, pol_type=None,
        temperature=1.0, scatter_axis=None,
        backend=None, backend_kws=None):
    """
    Calculate X-ray absorption spectra through a numerical backend.

    Backend-neutral physical arguments are passed directly. Numerical
    controls such as the SciPy Lanczos dimension belong in ``backend_kws``.
    When ``backend`` is omitted, it is inferred from ``hmat_n`` and the
    transition operators.
    """
    backend_name = (
        backend if backend is not None else _infer_backend(hmat_n, *trans_op)
    )
    match backend_name:
        case 'scipy':
            xas_backend = scipy_backend.xas_scipy
        case 'dense':
            xas_backend = scipy_backend.xas_dense
        case 'petsc':
            xas_backend = petsc_backend.xas_petsc
        case _:
            raise ValueError(
                "Unknown backend {!r}; expected 'scipy', 'dense', or "
                "'petsc'".format(backend_name)
            )

    return xas_backend(
        eval_i, evec_i, hmat_n, trans_op, ominc,
        gamma_c=gamma_c,
        thin=thin,
        phi=phi,
        pol_type=pol_type,
        temperature=temperature,
        scatter_axis=scatter_axis,
        backend_kws=backend_kws,
    )


def rixs(eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss, *,
         gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
         pol_type=None, temperature=1.0, scatter_axis=None,
         skip_gs=False, return_poles=False, backend=None, backend_kws=None):
    """
    Calculate resonant inelastic X-ray scattering spectra through a backend.

    Backend-neutral physical arguments are passed directly. Numerical
    controls such as SciPy Lanczos and GMRES settings belong in
    ``backend_kws``. When ``backend`` is omitted, it is inferred from the
    Hamiltonians and transition operators.

    If ``skip_gs`` is true, transitions into the retained initial-state
    subspace are omitted from the final-state spectrum, matching the legacy
    ``skip_gs`` behavior.
    """
    backend_name = (
        backend
        if backend is not None
        else _infer_backend(hmat_i, hmat_n, *trans_op)
    )
    match backend_name:
        case 'scipy':
            rixs_backend = scipy_backend.rixs_scipy
        case 'dense':
            rixs_backend = scipy_backend.rixs_dense
        case 'petsc':
            rixs_backend = petsc_backend.rixs_petsc
        case _:
            raise ValueError(
                "Unknown backend {!r}; expected 'scipy', 'dense', or "
                "'petsc'".format(backend_name)
            )

    return rixs_backend(
        eval_i, evec_i, hmat_i, hmat_n, trans_op, ominc, eloss,
        gamma_c=gamma_c,
        gamma_f=gamma_f,
        thin=thin,
        thout=thout,
        phi=phi,
        pol_type=pol_type,
        temperature=temperature,
        scatter_axis=scatter_axis,
        skip_gs=skip_gs,
        return_poles=return_poles,
        backend_kws=backend_kws,
    )


def _warn_legacy(name, replacement):
    """Emit the standard warning for a deprecated legacy solver workflow."""
    warnings.warn(
        "{} is deprecated; use {} instead.".format(name, replacement),
        DeprecationWarning,
        stacklevel=3,
    )


# -----------------------------------------------------------------------------
# Legacy dense-Python and Fortran workflows
# -----------------------------------------------------------------------------


def ed_1v1c_py(shell_name, *, shell_level=None, v_soc=None, c_soc=0,
               v_noccu=1, slater=None, ext_B=None, on_which='spin',
               v_cfmat=None, v_othermat=None, loc_axis=None, verbose=0):
    """
    Perform ED for the case of two atomic shells, one valence plus one Core
    shell with pure Python solver.
    For example, for Ni-:math:`L_3` edge RIXS, they are 3d valence and 2p core shells.

    It will use scipy.linalag.eigh to exactly diagonalize both the initial and intermediate
    Hamiltonians to get all the eigenvalues and eigenvectors, and the transition operators
    will be built in the many-body eigenvector basis.

    This solver is only suitable for small size of Hamiltonian, typically the dimension
    of both initial and intermediate Hamiltonian are smaller than 10,000.

    Parameters
    ----------
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') indicates a :math:`L_3` edge transition from
        core :math:`p_{3/2}` shell to valence :math:`d` shell.
    shell_level: tuple of two float numbers
        Energy level of valence (1st element) and core (2nd element) shells.

        They will be set to zero if not provided.
    v_soc: tuple of two float numbers
        Spin-orbit coupling strength of valence electrons, for the initial (1st element)
        and intermediate (2nd element) Hamiltonians.

        They will be set to zero if not provided.
    c_soc: a float number
        Spin-orbit coupling strength of core electrons.
    v_noccu: int number
        Number of electrons in valence shell.
    slater: tuple of two lists
        Slater integrals for initial (1st list) and intermediate (2nd list) Hamiltonians.
        The order of the elements in each list should be like this:

        [FX_vv, FX_vc, GX_vc, FX_cc],

        where X are integers with ascending order, it can be X=0, 2, 4, 6 or X=1, 3, 5.
        One can ignore all the continuous zeros at the end of the list.

        For example, if the full list is: [F0_dd, F2_dd, F4_dd, 0, F2_dp, 0, 0, 0, 0], one can
        just provide [F0_dd, F2_dd, F4_dd, 0, F2_dp]

        All the Slater integrals will be set to zero if slater=None.
    ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis.

        They will be set to zero if not provided.
    on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
    v_cfmat: 2d complex array
        Crystal field splitting Hamiltonian of valence electrons. The dimension and the orbital
        order should be consistent with the type of valence shell.

        They will be zeros if not provided.
    v_othermat: 2d complex array
        Other possible Hamiltonian of valence electrons. The dimension and the orbital order
        should be consistent with the type of valence shell.

        They will be zeros if not provided.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    verbose: int
        Level of writting data to files. Hopping matrices, Coulomb tensors, eigvenvalues
        will be written if verbose > 0.

    Returns
    -------
    eval_i:  1d float array
        The eigenvalues of initial Hamiltonian.
    eval_n: 1d float array
        The eigenvalues of intermediate Hamiltonian.
    trans_op: 3d complex array
        The matrices of transition operators in the eigenvector basis.
        Their components are defined with respect to the global :math:`xyz`-axis.
    """
    _warn_legacy(
        'ed_1v1c_py',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed',
    )
    print("edrixs >>> Running ED ...")
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    info_shell = info_atomic_shell()

    # Quantum numbers of angular momentum
    v_orbl = info_shell[v_name][0]

    # number of orbitals including spin degree of freedom
    v_norb = info_shell[v_name][1]
    c_norb = info_shell[c_name][1]

    # total number of orbitals
    ntot = v_norb + c_norb

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)

    # Coulomb interaction
    # Get the names of all the required slater integrals
    slater_name = slater_integrals_name((v_name, c_name), ('v', 'c'))
    nslat = len(slater_name)

    slater_i = np.zeros(nslat, dtype=float)
    slater_n = np.zeros(nslat, dtype=float)

    if slater is not None:
        if nslat > len(slater[0]):
            slater_i[0:len(slater[0])] = slater[0]
        else:
            slater_i[:] = slater[0][0:nslat]
        if nslat > len(slater[1]):
            slater_n[0:len(slater[1])] = slater[1]
        else:
            slater_n[:] = slater[1][0:nslat]

    # print summary of slater integrals
    print()
    print("    Summary of Slater integrals:")
    print("    ------------------------------")
    print("    Terms,   Initial Hamiltonian,  Intermediate Hamiltonian")
    for i in range(nslat):
        print("    ", slater_name[i], ":  {:20.10f}{:20.10f}".format(slater_i[i], slater_n[i]))
    print()

    case = v_name + c_name
    umat_i = get_umat_slater(case, *slater_i)
    umat_n = get_umat_slater(case, *slater_n)

    if verbose > 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

    # SOC
    if v_soc is not None:
        emat_i[0:v_norb, 0:v_norb] += atom_hsoc(v_name, v_soc[0])
        emat_n[0:v_norb, 0:v_norb] += atom_hsoc(v_name, v_soc[1])

    # when the core-shell is any of p12, p32, d32, d52, f52, f72,
    # do not need to add SOC for core shell
    if c_name in ['p', 'd', 'f']:
        emat_n[v_norb:ntot, v_norb:ntot] += atom_hsoc(c_name, c_soc)

    # crystal field
    if v_cfmat is not None:
        emat_i[0:v_norb, 0:v_norb] += np.array(v_cfmat)
        emat_n[0:v_norb, 0:v_norb] += np.array(v_cfmat)

    # other hopping matrix
    if v_othermat is not None:
        emat_i[0:v_norb, 0:v_norb] += np.array(v_othermat)
        emat_n[0:v_norb, 0:v_norb] += np.array(v_othermat)

    # energy of shells
    if shell_level is not None:
        emat_i[0:v_norb, 0:v_norb] += np.eye(v_norb) * shell_level[0]
        emat_i[v_norb:ntot, v_norb:ntot] += np.eye(c_norb) * shell_level[1]
        emat_n[0:v_norb, 0:v_norb] += np.eye(v_norb) * shell_level[0]
        emat_n[v_norb:ntot, v_norb:ntot] += np.eye(c_norb) * shell_level[1]

    # external magnetic field
    if v_name == 't2g':
        lx, ly, lz = get_lx(1, True), get_ly(1, True), get_lz(1, True)
        sx, sy, sz = get_sx(1), get_sy(1), get_sz(1)
        lx, ly, lz = -lx, -ly, -lz
    else:
        lx, ly, lz = get_lx(v_orbl, True), get_ly(v_orbl, True), get_lz(v_orbl, True)
        sx, sy, sz = get_sx(v_orbl), get_sy(v_orbl), get_sz(v_orbl)

    if ext_B is not None:
        if on_which.strip() == 'spin':
            zeeman = ext_B[0] * (2 * sx) + ext_B[1] * (2 * sy) + ext_B[2] * (2 * sz)
        elif on_which.strip() == 'orbital':
            zeeman = ext_B[0] * lx + ext_B[1] * ly + ext_B[2] * lz
        elif on_which.strip() == 'both':
            zeeman = ext_B[0] * (lx + 2 * sx) + ext_B[1] * (ly + 2 * sy) + ext_B[2] * (lz + 2 * sz)
        else:
            raise Exception("Unknown value of on_which", on_which)
        emat_i[0:v_norb, 0:v_norb] += zeeman
        emat_n[0:v_norb, 0:v_norb] += zeeman

    if verbose > 0:
        write_emat(emat_i, 'hopping_i.in')
        write_emat(emat_n, 'hopping_n.in')

    basis_i = get_fock_bin_by_N(v_norb, v_noccu, c_norb, c_norb)
    basis_n = get_fock_bin_by_N(v_norb, v_noccu+1, c_norb, c_norb - 1)
    ncfg_i, ncfg_n = len(basis_i), len(basis_n)
    print("edrixs >>> Dimension of the initial Hamiltonian: ", ncfg_i)
    print("edrixs >>> Dimension of the intermediate Hamiltonian: ", ncfg_n)

    # Build many-body Hamiltonian in Fock basis
    print("edrixs >>> Building Many-body Hamiltonians ...")
    hmat_i = np.zeros((ncfg_i, ncfg_i), dtype=complex)
    hmat_n = np.zeros((ncfg_n, ncfg_n), dtype=complex)
    hmat_i[:, :] += two_fermion(emat_i, basis_i, basis_i)
    hmat_i[:, :] += four_fermion(umat_i, basis_i)
    hmat_n[:, :] += two_fermion(emat_n, basis_n, basis_n)
    hmat_n[:, :] += four_fermion(umat_n, basis_n)
    print("edrixs >>> Done !")

    # Do exact-diagonalization to get eigenvalues and eigenvectors
    print("edrixs >>> Exact Diagonalization of Hamiltonians ...")
    eval_i, evec_i = scipy.linalg.eigh(hmat_i)
    eval_n, evec_n = scipy.linalg.eigh(hmat_n)
    print("edrixs >>> Done !")

    if verbose > 0:
        write_tensor(eval_i, 'eval_i.dat')
        write_tensor(eval_n, 'eval_n.dat')

    # Build dipolar transition operators in local-xyz axis
    if loc_axis is not None:
        local_axis = np.array(loc_axis)
    else:
        local_axis = np.eye(3)
    tmp = get_trans_oper(case)
    npol, n, m = tmp.shape
    tmp_g = np.zeros((npol, n, m), dtype=complex)
    # Transform the transition operators to global-xyz axis
    # dipolar transition
    if npol == 3:
        for i in range(3):
            for j in range(3):
                tmp_g[i] += local_axis[i, j] * tmp[j]

    # quadrupolar transition
    elif npol == 5:
        alpha, beta, gamma = rmat_to_euler(local_axis)
        wignerD = get_wigner_dmat(4, alpha, beta, gamma)
        rotmat = np.dot(np.dot(tmat_r2c('d'), wignerD), np.conj(np.transpose(tmat_r2c('d'))))
        for i in range(5):
            for j in range(5):
                tmp_g[i] += rotmat[i, j] * tmp[j]
    else:
        raise Exception("Have NOT implemented this case: ", npol)

    tmp2 = np.zeros((npol, ntot, ntot), dtype=complex)
    trans_op = np.zeros((npol, ncfg_n, ncfg_i), dtype=complex)
    for i in range(npol):
        tmp2[i, 0:v_norb, v_norb:ntot] = tmp_g[i]
        trans_op[i] = two_fermion(tmp2[i], basis_n, basis_i)
        trans_op[i] = cb_op2(trans_op[i], evec_n, evec_i)

    print("edrixs >>> ED Done !")

    return eval_i, eval_n, trans_op


def xas_1v1c_py(eval_i, eval_n, trans_op, ominc, *, gamma_c=0.1, thin=1.0, phi=0,
                pol_type=None, gs_list=None, temperature=1.0, scatter_axis=None):

    """
    Calculate XAS for the case of one valence shell plus one core shell with Python solver.

    This solver is only suitable for small size of Hamiltonian, typically the dimension
    of both initial and intermediate Hamiltonian are smaller than 10,000.

    Parameters
    ----------
    eval_i: 1d float array
        The eigenvalues of the initial Hamiltonian.
    eval_n: 1d float array
        The eigenvalues of the intermediate Hamiltonian.
    trans_op: 3d complex array
        The transition operators in the eigenstates basis.
    ominc: 1d float array
        Incident energy of photon.
    gamma_c: a float number or a 1d float array with the same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    thin: float number
        The incident angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of the scattering axis: scatter_axis[:,0].
    pol_type: list of tuples
        Type of polarization, options can be:

        - ('linear', alpha), linear polarization, where alpha is the angle between the
          polarization vector and the scattering plane in radians.

        - ('left', 0), left circular polarization.

        - ('right', 0), right circular polarization.

        - ('isotropic', 0). isotropic polarization.

        It will set pol_type=[('isotropic', 0)] if not provided.
    gs_list: 1d list of ints
        The indices of initial states which will be used in XAS calculations.

        It will set gs_list=[0] if not provided.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    scatter_axis: 3*3 float array
        The local axis defining the scattering plane. The scattering plane is defined in
        the local :math:`zx`-plane.

        local :math:`x`-axis: scatter_axis[:,0]

        local :math:`y`-axis: scatter_axis[:,1]

        local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    xas: 2d float array
        The calculated XAS spectra. The 1st dimension is for the incident energy, and the
        2nd dimension is for different polarizations.
    """
    _warn_legacy(
        'xas_1v1c_py',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed/xas',
    )

    print("edrixs >>> Running XAS ...")
    n_om = len(ominc)
    npol, ncfg_n = trans_op.shape[0], trans_op.shape[1]
    if pol_type is None:
        pol_type = [('isotropic', 0)]
    if gs_list is None:
        gs_list = [0]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    xas = np.zeros((n_om, len(pol_type)), dtype=float)
    gamma_core = np.zeros(n_om, dtype=float)
    prob = boltz_dist([eval_i[i] for i in gs_list], temperature)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_om) * gamma_c
    else:
        gamma_core[:] = gamma_c

    kvec = unit_wavevector(thin, phi, scatter_axis, 'in')
    for i, om in enumerate(ominc):
        for it, (pt, alpha) in enumerate(pol_type):
            if pt.strip() not in ['left', 'right', 'linear', 'isotropic']:
                raise Exception("Unknown polarization type: ", pt)
            polvec = np.zeros(npol, dtype=complex)
            if pt.strip() == 'left' or pt.strip() == 'right' or pt.strip() == 'linear':
                pol = dipole_polvec_xas(thin, phi, alpha, scatter_axis, pt)
                if npol == 3:  # dipolar transition
                    polvec[:] = pol
                if npol == 5:  # quadrupolar transition
                    polvec[:] = quadrupole_polvec(pol, kvec)

            # loop over all the initial states
            for j, igs in enumerate(gs_list):
                if pt.strip() == 'isotropic':
                    for k in range(npol):
                        xas[i, it] += (
                            prob[j] * np.sum(np.abs(trans_op[k, :, igs])**2 * gamma_core[i] /
                                             np.pi / ((om - (eval_n[:] - eval_i[igs]))**2 +
                                             gamma_core[i]**2))
                        ) / npol
                else:
                    F_mag = np.zeros(ncfg_n, dtype=complex)
                    for k in range(npol):
                        F_mag += trans_op[k, :, igs] * polvec[k]
                    xas[i, it] += (
                        prob[j] * np.sum(np.abs(F_mag)**2 * gamma_core[i] / np.pi /
                                         ((om - (eval_n[:] - eval_i[igs]))**2 + gamma_core[i]**2))
                    )

    print("edrixs >>> XAS Done !")

    return xas


def rixs_1v1c_py(eval_i, eval_n, trans_op, ominc, eloss, *,
                 gamma_c=0.1, gamma_f=0.01, thin=1.0, thout=1.0, phi=0.0,
                 pol_type=None, gs_list=None, temperature=1.0, scatter_axis=None,
                 skip_gs=False):
    """
    Calculate RIXS for the case of one valence shell plus one core shell with Python solver.

    This solver is only suitable for small size of Hamiltonian, typically the dimension
    of both initial and intermediate Hamiltonian are smaller than 10,000.

    Parameters
    ----------
    eval_i: 1d float array
        The eigenvalues of the initial Hamiltonian.
    eval_n: 1d float array
        The eigenvalues of the intermediate Hamiltonian.
    trans_op: 3d complex array
        The transition operators in the eigenstates basis.
    ominc: 1d float array
        Incident energy of photon.
    eloss: 1d float array
        Energy loss.
    gamma_c: a float number or a 1d float array with same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    gamma_f: a float number or a 1d float array with same shape as eloss.
        The final states life-time broadening factor. It can be a constant value
        or energy loss dependent.
    thin: float number
        The incident angle of photon (in radian).
    thout: float number
        The scattered angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of scattering axis: scatter_axis[:,0].
    pol_type: list of 4-elements-tuples
        Type of polarizations. It has the following form:

        (str1, alpha, str2, beta)

        where, str1 (str2) can be 'linear', 'left', 'right', 'isotropic' and alpha (beta) is
        the angle (in radians) between the linear polarization vector and the scattering plane.

        If str1 (or str2) is 'isotropic' then the polarization vector projects equally
        along each axis and the other variables are ignored.

        It will set pol_type=[('linear', 0, 'linear', 0)] if not provided.
    gs_list: 1d list of ints
        The indices of initial states which will be used in RIXS calculations.

        It will set gs_list=[0] if not provided.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    scatter_axis: 3*3 float array
        The local axis defining the scattering plane. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be an identity matrix if not provided.
    skip_gs: bool
        If True, transitions to the ground state(s) (forming the elastic peak)
        are omitted from the calculation.

    Returns
    -------
    rixs: 3d float array
        The calculated RIXS spectra. The 1st dimension is for the incident energy,
        the 2nd dimension is for the energy loss and the 3rd dimension is for
        different polarizations.
    """
    _warn_legacy(
        'rixs_1v1c_py',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed/rixs',
    )

    print("edrixs >>> Running RIXS ... ")
    n_ominc = len(ominc)
    n_eloss = len(eloss)
    gamma_core = np.zeros(n_ominc, dtype=float)
    gamma_final = np.zeros(n_eloss, dtype=float)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_ominc) * gamma_c
    else:
        gamma_core[:] = gamma_c

    if np.isscalar(gamma_f):
        gamma_final[:] = np.ones(n_eloss) * gamma_f
    else:
        gamma_final[:] = gamma_f

    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if gs_list is None:
        gs_list = [0]
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    prob = boltz_dist([eval_i[i] for i in gs_list], temperature)
    rixs = np.zeros((len(ominc), len(eloss), len(pol_type)), dtype=float)
    npol, n, m = trans_op.shape
    trans_emi = np.zeros((npol, m, n), dtype=np.complex128)
    for i in range(npol):
        trans_emi[i] = np.conj(np.transpose(trans_op[i]))
    polvec_i = np.zeros(npol, dtype=complex)
    polvec_f = np.zeros(npol, dtype=complex)

    # Calculate RIXS
    for i, om in enumerate(ominc):
        F_fi = scattering_mat(eval_i, eval_n, trans_op[:, :, 0:max(gs_list)+1],
                              trans_emi, om, gamma_core[i])

        for j, (it, alpha, jt, beta) in enumerate(pol_type):
            ei, ef = dipole_polvec_rixs(thin, thout, phi, alpha, beta,
                                        scatter_axis, (it, jt))
            if it.lower() == 'isotropic':
                ei = np.ones(3)/np.sqrt(3)                        # Powder spectrum
            if jt.lower() == 'isotropic':
                ef = np.ones(3)/np.sqrt(3)
            # dipolar transition
            if npol == 3:
                polvec_i[:] = ei
                polvec_f[:] = ef
            # quadrupolar transition
            elif npol == 5:
                ki = unit_wavevector(thin, phi, scatter_axis, direction='in')
                kf = unit_wavevector(thout, phi, scatter_axis, direction='out')
                polvec_i[:] = quadrupole_polvec(ei, ki)
                polvec_f[:] = quadrupole_polvec(ef, kf)
            else:
                raise Exception("Have NOT implemented this type of transition operators")
            # scattering magnitude with polarization vectors
            F_mag = np.zeros((len(eval_i), len(gs_list)), dtype=complex)
            for m in range(npol):
                for n in range(npol):
                    F_mag[:, :] += np.conj(polvec_f[m]) * F_fi[m, n] * polvec_i[n]

            fs_list = np.arange(len(eval_i))
            if skip_gs:
                fs_list = np.delete(fs_list, gs_list)
            for m, igs in enumerate(gs_list):
                for n in fs_list:
                    rixs[i, :, j] += (
                        prob[m] * np.abs(F_mag[n, igs])**2 * gamma_final / np.pi /
                        ((eloss - (eval_i[n] - eval_i[igs]))**2 + gamma_final**2)
                    )
    print("edrixs >>> RIXS Done !")

    return rixs


def ed_1v1c_fort(comm, shell_name, *, shell_level=None,
                 v_soc=None, c_soc=0, v_noccu=1, slater=None,
                 ext_B=None, on_which='spin',
                 v_cfmat=None, v_othermat=None,
                 do_ed=True, ed_solver=2, neval=1, nvector=1, ncv=3,
                 idump=False, maxiter=500, eigval_tol=1e-8, min_ndim=1000):
    """
    Perform ED for the case of one valence shell plus one Core-shell with Fortran ED solver.

    The hopping and Coulomb terms of both the initial and intermediate Hamiltonians will be
    constructed and written to files (hopping_i.in, hopping_n.in, coulomb_i.in and coulomb_n.in).
    Fock basis for the initial Hamiltonian will be written to file (fock_i.in).

    ED will be only performed on the initial Hamiltonian to find a few lowest eigenstates
    do_ed=True. Only input files will be written if do_ed=False.
    Due to large Hilbert space, the ed_fsolver written in Fortran will be called.
    mpi4py and a MPI environment (mpich or openmpi) are required to launch ed_fsolver.

    If do_ed=True, it will output the eigenvalues in file (eigvals.dat) and eigenvectors in files
    (eigvec.n), where n means the n-th eigenvectors. The eigvec.n files will be used later
    as the inputs for XAS and RIXS calculations.

    Parameters
    ----------
    comm: MPI_comm
        The MPI communicator from mpi4py.
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st strings can only be 's', 'p', 't2g', 'd', 'f'

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') may indicate a :math:`L_3` edge transition from
        core :math:`2p_{3/2}` shell to valence :math:`3d` shell for Ni.
    shell_level: tuple of two float numbers
        Energy level of valence (1st element) and core (2nd element) shells.

        They will be set to zero if not provided.
    v_soc: tuple of two float numbers
        Spin-orbit coupling strength of the valence shell,
        v1_soc[0] for the initial Hamiltonian, and
        v1_soc[1] for the intermediate Hamiltonian.

        They will be set to zero if not provided.
    c_soc: float number
        Spin-orbit coupling strength of core electrons.
    v_noccu: int number
        Total number of electrons in valence shells.
    slater: tuple of two lists
        Slater integrals for initial (1st list) and intermediate (2nd list) Hamiltonians.
        The order of the elements in each list should be like this:

        [FX_vv, FX_vc, GX_vc, FX_cc],

        where X are integers with ascending order, it can be X=0, 2, 4, 6 or X=1, 3, 5.
        One can ignore all the continuous zeros at the end of the list.

        For example, if the full list is: [F0_dd, F2_dd, F4_dd, 0, F2_dp, 0, 0, 0, 0], one can
        just provide [F0_dd, F2_dd, F4_dd, 0, F2_dp]

        All the Slater integrals will be set to zero if slater==None.
    ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis applied
        on the valence shell.

        It will be set to zeros if not provided.
    on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
    v_cfmat: 2d complex array
        Crystal field splitting Hamiltonian of the valence shell. The dimension and the orbital
        order should be consistent with the type of the valence shell.

        They will be zeros if not provided.
    v_othermat: 2d complex array
        Other possible Hamiltonian of the valence shell. The dimension and the orbital order
        should be consistent with the type of the valence shell.

        They will be zeros if not provided.
    do_ed: logical
        If do_end=True, diagonalize the Hamitlonian to find a few lowest eigenstates, return the
        eigenvalues and density matirx, and write the eigenvectors in files eigvec.n, otherwise,
        just write out the input files, do not perform the ED.
    ed_solver: int
        Type of ED solver, options can be 0, 1, 2

        - 0: use Lapack to fully diagonalize Hamiltonian to get all the eigenvalues.

        - 1: use standard Lanczos algorithm to find only a few lowest eigenvalues,
          no re-orthogonalization has been applied, so it is not very accurate.

        - 2: use parallel version of Arpack library to find a few lowest eigenvalues,
          it is accurate and is the recommeded choice in real calculations of XAS and RIXS.
    neval: int
        Number of eigenvalues to be found. For ed_solver=2, the value should not be too small,
        neval > 10 is usually a safe value.
    nvector: int
        Number of eigenvectors to be found and written into files.
    ncv: int
        Used for ed_solver=2, it should be at least ncv > neval + 2. Usually, set it a little
        bit larger than neval, for example, set ncv=200 when neval=100.
    idump: logical
        Whether to dump the eigenvectors to files "eigvec.n", where n means the n-th vectors.
    maxiter: int
        Maximum number of iterations in finding all the eigenvalues, used for ed_solver=1, 2.
    eigval_tol: float
        The convergence criteria of eigenvalues, used for ed_solver=1, 2.
    min_ndim: int
        The minimum dimension of the Hamiltonian when the ed_solver=1, 2 can be used, otherwise,
        ed_solver=1 will be used.

    Returns
    -------
    eval_i: 1d float array, shape=(neval, )
        The eigenvalues of initial Hamiltonian.
    denmat: 2d complex array, shape=(nvector, v_norb, v_norb))
        The density matrix in the eigenstates.
    """
    _warn_legacy(
        'ed_1v1c_fort',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    names = (v_name, 'empty', c_name)
    if shell_level is not None:
        levels = (shell_level[0], 0, shell_level[1])
    else:
        levels = None

    eval_i, denmat = _ed_1or2_valence_1core(
        comm, names, shell_level=levels, v1_soc=v_soc, c_soc=c_soc,
        v_tot_noccu=v_noccu, slater=slater, v1_ext_B=ext_B,
        v1_on_which=on_which, v1_cfmat=v_cfmat, v1_othermat=v_othermat,
        do_ed=do_ed, ed_solver=ed_solver, neval=neval, nvector=nvector,
        ncv=ncv, idump=idump, maxiter=maxiter, eigval_tol=eigval_tol,
        min_ndim=min_ndim
    )

    return eval_i, denmat


def xas_1v1c_fort(comm, shell_name, ominc, *, gamma_c=0.1,
                  v_noccu=1, thin=1.0, phi=0, pol_type=None,
                  num_gs=1, nkryl=200, temperature=1.0,
                  loc_axis=None, scatter_axis=None):
    """
    Calculate XAS for the case with one valence shells plus one core shell with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') may indicate a :math:`L_3` edge transition from
        core :math:`2p_{3/2}` shell to valence :math:`3d` shell for Ni.
    ominc: 1d float array
        Incident energy of photon.
    gamma_c: a float number or a 1d float array with the same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    v_noccu: int
        Total occupancy of valence shells.
    thin: float number
        The incident angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of the local scattering axis: scatter_axis[:,0].
    pol_type: list of tuples
        Type of polarization, options can be:

        - ('linear', alpha), linear polarization, where alpha is the angle between the
          polarization vector and the scattering plane in radians.

        - ('left', 0), left circular polarization.

        - ('right', 0), right circular polarization.

        - ('isotropic', 0). isotropic polarization.

        It will set pol_type=[('isotropic', 0)] if not provided.
    num_gs: int
        Number of initial states used in XAS calculations.
    nkryl: int
        Maximum number of poles obtained.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    xas: 2d array, shape=(len(ominc), len(pol_type))
        The calculated XAS spectra. The first dimension is for ominc, and the second dimension
        if for different polarizations.
    poles: list of dict, shape=(len(pol_type), )
        The calculated XAS poles for different polarizations.
    """
    _warn_legacy(
        'xas_1v1c_fort',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed/xas',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']

    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    names = (v_name, 'empty', c_name)

    xas, poles = _xas_1or2_valence_1core(
        comm, names, ominc, gamma_c=gamma_c, v_tot_noccu=v_noccu,
        trans_to_which=1, thin=thin, phi=phi, pol_type=pol_type,
        num_gs=num_gs, nkryl=nkryl, temperature=temperature,
        loc_axis=loc_axis, scatter_axis=scatter_axis
    )

    return xas, poles


def rixs_1v1c_fort(comm, shell_name, ominc, eloss, *, gamma_c=0.1, gamma_f=0.1,
                   v_noccu=1, thin=1.0, thout=1.0, phi=0, pol_type=None,
                   num_gs=1, nkryl=200, linsys_max=500, linsys_tol=1e-8,
                   temperature=1.0, loc_axis=None, scatter_axis=None):
    """
    Calculate RIXS for the case with one valence shell plus one core shell with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') may indicate a :math:`L_3` edge transition from
        core :math:`2p_{3/2}` shell to valence :math:`3d` shell for Ni.
    ominc: 1d float array
        Incident energy of photon.
    eloss: 1d float array
        Energy loss.
    gamma_c: a float number or a 1d float array with same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    gamma_f: a float number or a 1d float array with same shape as eloss.
        The final states life-time broadening factor. It can be a constant value
        or energy loss dependent.
    v_noccu: int
        Total occupancy of valence shells.
    thin: float number
        The incident angle of photon (in radian).
    thout: float number
        The scattered angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of scattering axis: scatter_axis[:,0].
    pol_type: list of 4-elements-tuples
        Type of polarizations. It has the following form:

        (str1, alpha, str2, beta)

        where, str1 (str2) can be 'linear', 'left', 'right', and alpha (beta) is
        the angle (in radians) between the linear polarization vector and the scattering plane.

        It will set pol_type=[('linear', 0, 'linear', 0)] if not provided.
    num_gs: int
        Number of initial states used in RIXS calculations.
    nkryl: int
        Maximum number of poles obtained.
    linsys_max: int
        Maximum iterations of solving linear equations.
    linsys_tol: float
        Convergence for solving linear equations.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    rixs: 3d float array, shape=(len(ominc), len(eloss), len(pol_type))
        The calculated RIXS spectra. The 1st dimension is for the incident energy,
        the 2nd dimension is for the energy loss and the 3rd dimension is for
        different polarizations.
    poles: 2d list of dict, shape=(len(ominc), len(pol_type))
        The calculated RIXS poles. The 1st dimension is for incident energy, and the
        2nd dimension is for different polarizations.
    """
    _warn_legacy(
        'rixs_1v1c_fort',
        'edrixs.models.model_1v1c + edrixs.solvers.get_ops/ed/rixs',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    names = (v_name, 'empty', c_name)
    rixs, poles = _rixs_1or2_valence_1core(
        comm, names, ominc, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        v_tot_noccu=v_noccu, trans_to_which=1, thin=thin,
        thout=thout, phi=phi, pol_type=pol_type, num_gs=num_gs, nkryl=nkryl,
        linsys_max=linsys_max, linsys_tol=linsys_tol, temperature=temperature,
        loc_axis=loc_axis, scatter_axis=loc_axis
    )

    return rixs, poles


def ed_2v1c_fort(comm, shell_name, *, shell_level=None,
                 v1_soc=None, v2_soc=None, c_soc=0, v_tot_noccu=1, slater=None,
                 v1_ext_B=None, v2_ext_B=None, v1_on_which='spin', v2_on_which='spin',
                 v1_cfmat=None, v2_cfmat=None, v1_othermat=None, v2_othermat=None,
                 hopping_v1v2=None, do_ed=True, ed_solver=2, neval=1, nvector=1, ncv=3,
                 idump=False, maxiter=500, eigval_tol=1e-8, min_ndim=1000):
    """
    Perform ED for the case of two valence shell plus one core-shell with Fortran solver.
    For example, for Ni :math:`K`-edge RIXS, :math:`1s\\rightarrow 4p` transition,
    the valence shells involved in RIXS are :math:`3d` and :math:`4p`.

    The hopping and Coulomb terms of both the initial and intermediate Hamiltonians will be
    constructed and written to files (hopping_i.in, hopping_n.in, coulomb_i.in and coulomb_n.in).
    Fock basis for the initial Hamiltonian will be written to file (fock_i.in).

    ED will be only performed on the initial Hamiltonian to find a few lowest eigenstates
    do_ed=True. Only input files will be written if do_ed=False.
    Due to large Hilbert space, the ed_fsolver written in Fortran will be called.
    mpi4py and a MPI environment (mpich or openmpi) are required to launch ed_fsolver.

    If do_ed=True, it will output the eigenvalues in file (eigvals.dat) and eigenvectors in files
    (eigvec.n), where n means the n-th eigenvectors. The eigvec.n files will be used later
    as the inputs for XAS and RIXS calculations.

    Parameters
    ----------
    comm: MPI_comm
        The MPI communicator from mpi4py.
    shell_name: tuple of three strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        1st (2nd) valence shell, and the 3rd one is for the core shell.

        - The 1st and 2nd strings can only be 's', 'p', 't2g', 'd', 'f'

        - The 3nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p', 's') may indicate a :math:`K` edge transition from
        core :math:`1s` shell to valence :math:`3d` and :math:`4p` shell for Ni.
    shell_level: tuple of three float numbers
        Energy level of valence (1st and 2nd elements) and core (3nd element) shells.

        They will be set to zero if not provided.
    v1_soc: tuple of two float numbers
        Spin-orbit coupling strength of the 1st valence shell,
        v1_soc[0] for the initial Hamiltonian, and
        v1_soc[1] for the intermediate Hamiltonian.

        They will be set to zero if not provided.
    v2_soc: tuple of two float numbers
        Spin-orbit coupling strength of the 2nd valence shell,
        v2_soc[0] for the initial Hamiltonian, and
        v2_soc[1] for the intermediate Hamiltonian.

        They will be set to zero if not provided.
    c_soc: float number
        Spin-orbit coupling strength of core electrons.
    v_tot_noccu: int number
        Total number of electrons in valence shells.
    slater: tuple of two lists
        Slater integrals for initial (1st list) and intermediate (2nd list) Hamiltonians.
        The order of the elements in each list should be like this:

        [FX_v1v1, FX_v1v2, GX_v1v2, FX_v2v2, FX_v1c, GX_v1c, FX_v2c, GX_v2c],

        where X are integers with ascending order, it can be X=0, 2, 4, 6 or X=1, 3, 5.
        One can ignore all the continuous zeros at the end of the list.

        For example, if the full list is: [F0_dd, F2_dd, F4_dd, 0, F2_dp, 0, 0, 0, 0], one can
        just provide [F0_dd, F2_dd, F4_dd, 0, F2_dp]

        All the Slater integrals will be set to zero if slater==None.
    v1_ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis applied
        on the 1st valence shell.

        It will be set to zeros if not provided.
    v2_ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis applied
        on the 2nd valence shell.

        It will be set to zeros if not provided.
    v1_on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
        For the 1st valence shell.
    v2_on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
        For the 2nd valence shell.
    v1_cfmat: 2d complex array
        Crystal field splitting Hamiltonian of the 1st valence shell. The dimension and the orbital
        order should be consistent with the type of the 1st valence shell.

        They will be zeros if not provided.
    v2_cfmat: 2d complex array
        Crystal field splitting Hamiltonian of the 2nd valence shell. The dimension and the orbital
        order should be consistent with the type of the 2nd valence shell.

        They will be zeros if not provided.
    v1_othermat: 2d complex array
        Other possible Hamiltonian of the 1st valence shell. The dimension and the orbital order
        should be consistent with the type of the 1st valence shell.

        They will be zeros if not provided.
    v2_othermat: 2d complex array
        Other possible Hamiltonian of the 2nd valence shell. The dimension and the orbital order
        should be consistent with the type of the 2nd valence shell.

        They will be zeros if not provided.
    hopping_v1v2: 2d complex array
        Hopping between the two valence shells. The 1st-index (2nd-index) is the 1st (2nd)
        valence shell.

        They will be zeros if not provided.
    do_ed: logical
        If do_end=True, diagonalize the Hamitlonian to find a few lowest eigenstates, return the
        eigenvalues and density matirx, and write the eigenvectors in files eigvec.n, otherwise,
        just write out the input files, do not perform the ED.
    ed_solver: int
        Type of ED solver, options can be 0, 1, 2

        - 0: use Lapack to fully diagonalize Hamiltonian to get all the eigenvalues.

        - 1: use standard Lanczos algorithm to find only a few lowest eigenvalues,
          no re-orthogonalization has been applied, so it is not very accurate.

        - 2: use parallel version of Arpack library to find a few lowest eigenvalues,
          it is accurate and is the recommeded choice in real calculations of XAS and RIXS.
    neval: int
        Number of eigenvalues to be found. For ed_solver=2, the value should not be too small,
        neval > 10 is usually a safe value.
    nvector: int
        Number of eigenvectors to be found and written into files.
    ncv: int
        Used for ed_solver=2, it should be at least ncv > neval + 2. Usually, set it a little
        bit larger than neval, for example, set ncv=200 when neval=100.
    idump: logical
        Whether to dump the eigenvectors to files "eigvec.n", where n means the n-th vectors.
    maxiter: int
        Maximum number of iterations in finding all the eigenvalues, used for ed_solver=1, 2.
    eigval_tol: float
        The convergence criteria of eigenvalues, used for ed_solver=1, 2.
    min_ndim: int
        The minimum dimension of the Hamiltonian when the ed_solver=1, 2 can be used, otherwise,
        ed_solver=1 will be used.

    Returns
    -------
    eval_i: 1d float array, shape=(neval, )
        The eigenvalues of initial Hamiltonian.
    denmat: 2d complex array, shape=(nvector, v1v2_norb, v1v2_norb))
        The density matrix in the eigenstates.
    """
    _warn_legacy(
        'ed_2v1c_fort',
        'edrixs.models.model_2v1c + edrixs.solvers.get_ops/ed',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()
    if v1_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v1_name)
    if v2_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v2_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    eval_i, denmat = _ed_1or2_valence_1core(
        comm, shell_name, shell_level=shell_level, v1_soc=v1_soc, v2_soc=v2_soc, c_soc=c_soc,
        v_tot_noccu=v_tot_noccu, slater=slater, v1_ext_B=v1_ext_B, v2_ext_B=v2_ext_B,
        v1_on_which=v1_on_which, v2_on_which=v2_on_which, v1_cfmat=v1_cfmat,
        v2_cfmat=v2_cfmat, v1_othermat=v1_othermat, v2_othermat=v2_othermat,
        hopping_v1v2=hopping_v1v2, do_ed=do_ed, ed_solver=ed_solver, neval=neval,
        nvector=nvector, ncv=ncv, idump=idump, maxiter=maxiter, eigval_tol=eigval_tol,
        min_ndim=min_ndim
    )
    return eval_i, denmat


def xas_2v1c_fort(comm, shell_name, ominc, *, gamma_c=0.1,
                  v_tot_noccu=1, trans_to_which=1, thin=1.0, phi=0,
                  pol_type=None, num_gs=1, nkryl=200, temperature=1.0,
                  loc_axis=None, scatter_axis=None):
    """
    Calculate XAS for the case with two valence shells plus one core shell with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of three strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        1st (2nd) valence shell, and the 3rd one is for the core shell.

        - The 1st and 2nd strings can only be 's', 'p', 't2g', 'd', 'f',

        - The 3nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p', 's') may indicate a :math:`K` edge transition from
        core :math:`1s` shell to valence :math:`3d` and :math:`4p` shell for Ni.
    ominc: 1d float array
        Incident energy of photon.
    gamma_c: a float number or a 1d float array with the same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    v_tot_noccu: int
        Total occupancy of valence shells.
    trans_to_which: int
        Photon transition to which valence shell.

        - 1: to 1st valence shell,

        - 2: to 2nd valence shell.
    thin: float number
        The incident angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of the local scattering axis: scatter_axis[:,0].
    pol_type: list of tuples
        Type of polarization, options can be:

        - ('linear', alpha), linear polarization, where alpha is the angle between the
          polarization vector and the scattering plane in radians.

        - ('left', 0), left circular polarization.

        - ('right', 0), right circular polarization.

        - ('isotropic', 0). isotropic polarization.

        It will set pol_type=[('isotropic', 0)] if not provided.
    num_gs: int
        Number of initial states used in XAS calculations.
    nkryl: int
        Maximum number of poles obtained.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    xas: 2d array, shape=(len(ominc), len(pol_type))
        The calculated XAS spectra. The first dimension is for ominc, and the second dimension
        if for different polarizations.
    poles: list of dict, shape=(len(pol_type), )
        The calculated XAS poles for different polarizations.
    """
    _warn_legacy(
        'xas_2v1c_fort',
        'edrixs.models.model_2v1c + edrixs.solvers.get_ops/ed/xas',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']

    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()
    if v1_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v1_name)
    if v2_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v2_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    xas, poles = _xas_1or2_valence_1core(
        comm, shell_name, ominc, gamma_c=gamma_c, v_tot_noccu=v_tot_noccu,
        trans_to_which=trans_to_which, thin=thin, phi=phi, pol_type=pol_type,
        num_gs=num_gs, nkryl=nkryl, temperature=temperature,
        loc_axis=loc_axis, scatter_axis=scatter_axis
    )

    return xas, poles


def rixs_2v1c_fort(comm, shell_name, ominc, eloss, *, gamma_c=0.1, gamma_f=0.1,
                   v_tot_noccu=1, trans_to_which=1, thin=1.0, thout=1.0, phi=0,
                   pol_type=None, num_gs=1, nkryl=200, linsys_max=500, linsys_tol=1e-8,
                   temperature=1.0, loc_axis=None, scatter_axis=None):
    """
    Calculate RIXS for the case with 2 valence shells plus 1 core shell.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of three strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        1st (2nd) valence shell, and the 3rd one is for the core shell.

        - The 1st and 2nd strings can only be 's', 'p', 't2g', 'd', 'f',

        - The 3nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p', 's') may indicate a :math:`K` edge transition from
        core :math:`1s` shell to valence :math:`3d` and :math:`4p` shell for Ni.
    ominc: 1d float array
        Incident energy of photon.
    eloss: 1d float array
        Energy loss.
    gamma_c: a float number or a 1d float array with same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    gamma_f: a float number or a 1d float array with same shape as eloss.
        The final states life-time broadening factor. It can be a constant value
        or energy loss dependent.
    v_tot_noccu: int
        Total occupancy of valence shells.
    trans_to_which: int
        Photon transition to which valence shell.

        - 1: to 1st valence shell,

        - 2: to 2nd valence shell.
    thin: float number
        The incident angle of photon (in radian).
    thout: float number
        The scattered angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of scattering axis: scatter_axis[:,0].
    pol_type: list of 4-elements-tuples
        Type of polarizations. It has the following form:

        (str1, alpha, str2, beta)

        where, str1 (str2) can be 'linear', 'left', 'right', and alpha (beta) is
        the angle (in radians) between the linear polarization vector and the scattering plane.

        It will set pol_type=[('linear', 0, 'linear', 0)] if not provided.
    num_gs: int
        Number of initial states used in RIXS calculations.
    nkryl: int
        Maximum number of poles obtained.
    linsys_max: int
        Maximum iterations of solving linear equations.
    linsys_tol: float
        Convergence for solving linear equations.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    rixs: 3d float array, shape=(len(ominc), len(eloss), len(pol_type))
        The calculated RIXS spectra. The 1st dimension is for the incident energy,
        the 2nd dimension is for the energy loss and the 3rd dimension is for
        different polarizations.
    poles: 2d list of dict, shape=(len(ominc), len(pol_type))
        The calculated RIXS poles. The 1st dimension is for incident energy, and the
        2nd dimension is for different polarizations.
    """
    _warn_legacy(
        'rixs_2v1c_fort',
        'edrixs.models.model_2v1c + edrixs.solvers.get_ops/ed/rixs',
    )
    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()
    if v1_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v1_name)
    if v2_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v2_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    rixs, poles = _rixs_1or2_valence_1core(
        comm, shell_name, ominc, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        v_tot_noccu=v_tot_noccu, trans_to_which=trans_to_which, thin=thin,
        thout=thout, phi=phi, pol_type=pol_type, num_gs=num_gs, nkryl=nkryl,
        linsys_max=linsys_max, linsys_tol=linsys_tol, temperature=temperature,
        loc_axis=loc_axis, scatter_axis=loc_axis
    )

    return rixs, poles


def ed_siam_fort(comm, shell_name, nbath, *, siam_type=0, v_noccu=1, static_core_pot=0, c_level=0,
                 c_soc=0, trans_c2n=None, imp_mat=None, imp_mat_n=None, bath_level=None,
                 bath_level_n=None, hyb=None, hyb_n=None, hopping=None, hopping_n=None,
                 slater=None, ext_B=None, on_which='spin', do_ed=0, ed_solver=2, neval=1,
                 nvector=1, ncv=3, idump=False, maxiter=1000, eigval_tol=1e-8, min_ndim=1000):
    """
    Find the ground state of the initial Hamiltonian of a Single Impuirty Anderson Model (SIAM),
    and also prepare input files, *hopping_i.in*, *hopping_n.in*, *coulomb_i.in*, *coulomb_n.in*
    for following XAS and RIXS calculations.

    Parameters
    ----------
    comm: MPI_Comm
        MPI Communicator
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') indicates a :math:`L_3` edge transition from
        core :math:`p_{3/2}` shell to valence :math:`d` shell.
    nbath: int
        Number of bath sites.
    siam_type: int
        Type of SIAM Hamiltonian,

        - 0: diagonal hybridization function, parameterized by *imp_mat*, *bath_level* and *hyb*

        - 1: general hybridization function, parameterized by matrix *hopping*

        if *siam_type=0*, only *imp_mat*, *bath_level* and *hyb* are required,
        if *siam_type=1*, only *hopping* is required.
    v_noccu: int
        Number of total occupancy of impurity and baths orbitals, required when do_ed=1, 2
    static_core_pot: float
        Static core hole potential.
    c_level: float
        Energy level of core shell.
    c_soc: float
        Spin-orbit coupling strength of core electrons.
    trans_c2n: 2d complex array
        The transformation matrix from the spherical harmonics basis to the basis on which
        the `imp_mat` and hybridization function (`bath_level`, `hyb`, `hopping`) are defined.
    imp_mat: 2d complex array
        Impurity matrix for the impurity site, including CF or SOC, for siam_type=0
        and the initial configurations.
    imp_mat_n: 2d complex array
        Impurity matrix for the impurity site, including CF or SOC, for siam_type=0
        and the intermediate configurations. If imp_mat_n=None, imp_mat will be used.
    bath_level: 2d complex array
        Energy level of bath sites, 1st (2nd) dimension is for different bath sites (orbitals),
        for siam_type=0 and the initial configurations.
    bath_level_n: 2d complex array
        Energy level of bath sites, 1st (2nd) dimension is for different bath sites (orbitals),
        for siam_type=0 and the intermediate configurations. If bath_level_n=None,
        bath_level will be used.
    hyb: 2d complex array
        Hybridization strength of bath sites, 1st (2nd) dimension is for different bath
        sites (orbitals), for siam_type=0 and the initial configurations.
    hyb_n: 2d complex array
        Hybridization strength of bath sites, 1st (2nd) dimension is for different bath
        sites (orbitals), for siam_type=0 and the intermediate configurations.
        If hyb_n=None, hyb will be used.
    hopping: 2d complex array
        General hopping matrix when siam_type=1, including imp_mat and hybridization functions,
        for siam_type=1 and the initial configurations.
    hopping_n: 2d complex array
        General hopping matrix when siam_type=1, including imp_mat and hybridization functions,
        for siam_type=1 and the intermediate configurations. If hopping_n=None,
        hopping will be used.
    slater: tuple of two lists
        Slater integrals for initial (1st list) and intermediate (2nd list) Hamiltonians.
        The order of the elements in each list should be like this:

        [FX_vv, FX_vc, GX_vc, FX_cc],

        where X are integers with ascending order, it can be X=0, 2, 4, 6 or X=1, 3, 5.
        One can ignore all the continuous zeros at the end of the list.

        For example, if the full list is: [F0_dd, F2_dd, F4_dd, 0, F2_dp, 0, 0, 0, 0], one can
        just provide [F0_dd, F2_dd, F4_dd, 0, F2_dp]

        All the Slater integrals will be set to zero if slater=None.
    ext_B: tuple of three float numbers
        Vector of external magnetic field with respect to global :math:`xyz`-axis.

        They will be set to zero if not provided.
    on_which: string
        Apply Zeeman exchange field on which sector. Options are 'spin', 'orbital' or 'both'.
    do_ed: int
        - 0: First, search the ground state in different subspaces of total occupancy
          :math:`N` with ed_solver=1, and then do a more accurate ED in the subspace
          :math:`N` where the ground state lies to find a few lowest eigenstates, return
          the eigenvalues and density matirx, and write the eigenvectors in files eigvec.n

        - 1: Only do ED for given occupancy number *v_noccu*, return eigenvalues and
          density matrix, write eigenvectors to files eigvec.n

        - 2: Do not do ED, only write parameters into files: *hopping_i.in*, *hopping_n.in*,
          *coulomb_i.in*, *coulomb_n.in* for later XAS or RIXS calculations.
    ed_solver: int
        Type of ED solver, options can be 0, 1, 2

        - 0: use Lapack to fully diagonalize Hamiltonian to get all the eigenvalues.

        - 1: use standard Lanczos algorithm to find only a few lowest eigenvalues,
          no re-orthogonalization has been applied, so it is not very accurate.

        - 2: use parallel version of Arpack library to find a few lowest eigenvalues,
          it is accurate and is the recommeded choice in real calculations of XAS and RIXS.
    neval: int
        Number of eigenvalues to be found. For ed_solver=2, the value should not be too small,
        neval > 10 is usually a safe value.
    nvector: int
        Number of eigenvectors to be found and written into files.
    ncv: int
        Used for ed_solver=2, it should be at least ncv > neval + 2. Usually, set it a little
        bit larger than neval, for example, set ncv=200 when neval=100.
    idump: logical
        Whether to dump the eigenvectors to files "eigvec.n", where n means the n-th vectors.
    maxiter: int
        Maximum number of iterations in finding all the eigenvalues, used for ed_solver=1, 2.
    eigval_tol: float
        The convergence criteria of eigenvalues, used for ed_solver=1, 2.
    min_ndim: int
        The minimum dimension of the Hamiltonian when the ed_solver=1, 2 can be used, otherwise,
        ed_solver=1 will be used.

    Returns
    -------
    eval_i: 1d float array
        Eigenvalues of initial Hamiltonian.
    denmat: 2d complex array
        Density matrix.
    noccu_gs: int
        Occupancy of the ground state.
    """
    _warn_legacy(
        'ed_siam_fort',
        'edrixs.models.model_siam + edrixs.solvers.get_ops/ed',
    )
    from .fedrixs import ed_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()
    if rank == 0:
        print("edrixs >>> Running ED ...", flush=True)

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    info_shell = info_atomic_shell()
    v_orbl = info_shell[v_name][0]
    v_norb = info_shell[v_name][1]
    c_norb = info_shell[c_name][1]
    ntot_v = v_norb * (nbath + 1)
    ntot = ntot_v + c_norb

    slater_name = slater_integrals_name((v_name, c_name), ('v', 'c'))
    nslat = len(slater_name)
    slater_i = np.zeros(nslat, dtype=float)
    slater_n = np.zeros(nslat, dtype=float)

    if slater is not None:
        if nslat > len(slater[0]):
            slater_i[0:len(slater[0])] = slater[0]
        else:
            slater_i[:] = slater[0][0:nslat]
        if nslat > len(slater[1]):
            slater_n[0:len(slater[1])] = slater[1]
        else:
            slater_n[:] = slater[1][0:nslat]

    # print summary of slater integrals
    if rank == 0:
        print(flush=True)
        print("    Summary of Slater integrals:", flush=True)
        print("    ------------------------------", flush=True)
        print("    Terms,  Initial Hamiltonian,  Intermediate Hamiltonian", flush=True)
        for i in range(nslat):
            print(
                "    ", slater_name[i],
                ":  {:20.10f}{:20.10f}".format(slater_i[i], slater_n[i]), flush=True
            )
        print(flush=True)

    umat_tmp_i = get_umat_slater(v_name + c_name, *slater_i)
    umat_tmp_n = get_umat_slater(v_name + c_name, *slater_n)

    umat_i = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)
    umat_n = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)

    indx = list(range(0, v_norb)) + [ntot_v + i for i in range(0, c_norb)]
    for i in range(v_norb+c_norb):
        for j in range(v_norb+c_norb):
            for k in range(v_norb+c_norb):
                for m in range(v_norb+c_norb):
                    umat_i[indx[i], indx[j], indx[k], indx[m]] = umat_tmp_i[i, j, k, m]
                    umat_n[indx[i], indx[j], indx[k], indx[m]] = umat_tmp_n[i, j, k, m]
    if rank == 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)
    # General hybridization function, including off-diagonal terms
    if siam_type == 1:
        if hopping is not None:
            emat_i[0:ntot_v, 0:ntot_v] += hopping
        if hopping_n is not None:
            emat_n[0:ntot_v, 0:ntot_v] += hopping_n
        elif hopping is not None:
            emat_n[0:ntot_v, 0:ntot_v] += hopping
    # Diagonal hybridization function
    elif siam_type == 0:
        # matrix (CF or SOC) for impuirty site
        if imp_mat is not None:
            emat_i[0:v_norb, 0:v_norb] += imp_mat
        if imp_mat_n is not None:
            emat_n[0:v_norb, 0:v_norb] += imp_mat_n
        elif imp_mat is not None:
            emat_n[0:v_norb, 0:v_norb] += imp_mat
        # bath levels
        if bath_level is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx = (i + 1) * v_norb + j
                    emat_i[indx, indx] += bath_level[i, j]
        if bath_level_n is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx = (i + 1) * v_norb + j
                    emat_n[indx, indx] += bath_level_n[i, j]
        elif bath_level is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx = (i + 1) * v_norb + j
                    emat_n[indx, indx] += bath_level[i, j]
        if hyb is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx1, indx2 = j, (i + 1) * v_norb + j
                    emat_i[indx1, indx2] += hyb[i, j]
                    emat_i[indx2, indx1] += np.conj(hyb[i, j])
        if hyb_n is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx1, indx2 = j, (i + 1) * v_norb + j
                    emat_n[indx1, indx2] += hyb_n[i, j]
                    emat_n[indx2, indx1] += np.conj(hyb_n[i, j])
        elif hyb is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    indx1, indx2 = j, (i + 1) * v_norb + j
                    emat_n[indx1, indx2] += hyb[i, j]
                    emat_n[indx2, indx1] += np.conj(hyb[i, j])
    else:
        raise Exception("Unknown siam_type: ", siam_type)

    if c_name in ['p', 'd', 'f']:
        emat_n[ntot_v:ntot, ntot_v:ntot] += atom_hsoc(c_name, c_soc)

    # static core potential
    emat_n[0:v_norb, 0:v_norb] -= np.eye(v_norb) * static_core_pot

    if trans_c2n is None:
        trans_c2n = np.eye(v_norb, dtype=complex)
    else:
        trans_c2n = np.array(trans_c2n)

    tmat = np.eye(ntot, dtype=complex)
    for i in range(nbath+1):
        off = i * v_norb
        tmat[off:off+v_norb, off:off+v_norb] = np.conj(np.transpose(trans_c2n))
    emat_i[:, :] = cb_op(emat_i, tmat)
    emat_n[:, :] = cb_op(emat_n, tmat)

    # zeeman field
    if v_name == 't2g':
        lx, ly, lz = get_lx(1, True), get_ly(1, True), get_lz(1, True)
        sx, sy, sz = get_sx(1), get_sy(1), get_sz(1)
        lx, ly, lz = -lx, -ly, -lz
    else:
        lx, ly, lz = get_lx(v_orbl, True), get_ly(v_orbl, True), get_lz(v_orbl, True)
        sx, sy, sz = get_sx(v_orbl), get_sy(v_orbl), get_sz(v_orbl)

    if ext_B is not None:
        if on_which.strip() == 'spin':
            zeeman = ext_B[0] * (2 * sx) + ext_B[1] * (2 * sy) + ext_B[2] * (2 * sz)
        elif on_which.strip() == 'orbital':
            zeeman = ext_B[0] * lx + ext_B[1] * ly + ext_B[2] * lz
        elif on_which.strip() == 'both':
            zeeman = ext_B[0] * (lx + 2 * sx) + ext_B[1] * (ly + 2 * sy) + ext_B[2] * (lz + 2 * sz)
        else:
            raise Exception("Unknown value of on_which", on_which)
        emat_i[0:v_norb, 0:v_norb] += zeeman
        emat_n[0:v_norb, 0:v_norb] += zeeman

    # Perform ED if necessary
    if do_ed == 1 or do_ed == 2:
        eval_shift = c_level * c_norb / v_noccu
        emat_i[0:ntot_v, 0:ntot_v] += np.eye(ntot_v) * eval_shift
        emat_n[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level
        if rank == 0:
            write_emat(emat_i, 'hopping_i.in')
            write_emat(emat_n, 'hopping_n.in')
            write_config(
                ed_solver=ed_solver, num_val_orbs=ntot_v, neval=neval, nvector=nvector, ncv=ncv,
                idump=idump, maxiter=maxiter, min_ndim=min_ndim, eigval_tol=eigval_tol
            )
            write_fock_dec_by_N(ntot_v, v_noccu, "fock_i.in")
        if do_ed == 1:
            if rank == 0:
                print("edrixs >>> do_ed=1, perform ED at noccu: ", v_noccu, flush=True)
            comm.Barrier()
            ed_fsolver(fcomm, rank, size)
            comm.Barrier()
            data = np.loadtxt('eigvals.dat', ndmin=2)
            eval_i = np.zeros(neval, dtype=float)
            eval_i[0:neval] = data[0:neval, 1]
            data = np.loadtxt('denmat.dat', ndmin=2)
            tmp = (nvector, ntot_v, ntot_v)
            denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
            return eval_i, denmat, v_noccu
        else:
            if rank == 0:
                print("edrixs >>> do_ed=2, Do not perform ED, only write files", flush=True)
            return None, None, None

    # Find the ground states by total occupancy N
    elif do_ed == 0:
        if rank == 0:
            print("edrixs >>> do_ed=0, serach ground state by total occupancy N", flush=True)
            flog = open('search_gs.log', 'w')
        write_emat(emat_i, 'hopping_i.in')
        res = []
        num_electron = ntot_v // 2
        noccu_gs = num_electron
        if rank == 0:
            write_config(
                ed_solver=1, num_val_orbs=ntot_v, neval=1, nvector=1, idump=False,
                maxiter=maxiter, min_ndim=min_ndim, eigval_tol=eigval_tol
            )
            write_fock_dec_by_N(ntot_v, num_electron, "fock_i.in")
        comm.Barrier()
        ed_fsolver(fcomm, rank, size)
        comm.Barrier()
        data = np.loadtxt('eigvals.dat', ndmin=2)
        eval_gs = data[0, 1]
        data = np.loadtxt('denmat.dat', ndmin=2)
        tmp = (1, ntot_v, ntot_v)
        denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
        imp_occu = np.sum(denmat[0].diagonal()[0:v_norb]).real
        res.append((num_electron, eval_gs, imp_occu))
        if rank == 0:
            print(num_electron, eval_gs, imp_occu, file=flog, flush=True)

        nplus_list = [num_electron + i + 1 for i in range(ntot_v // 2)]
        nminus_list = [num_electron - i - 1 for i in range(ntot_v // 2)]
        nplus_direction = True
        nminus_direction = True
        for i in range(ntot_v // 2):
            if nplus_direction:
                num_electron = nplus_list[i]
                if rank == 0:
                    write_fock_dec_by_N(ntot_v, num_electron, "fock_i.in")
                comm.Barrier()
                ed_fsolver(fcomm, rank, size)
                comm.Barrier()
                data = np.loadtxt('eigvals.dat', ndmin=2)
                eigval = data[0, 1]
                if eigval > eval_gs:
                    nplus_direction = False
                else:
                    nminus_direction = False
                    eval_gs = eigval
                    noccu_gs = num_electron
                data = np.loadtxt('denmat.dat', ndmin=2)
                tmp = (1, ntot_v, ntot_v)
                denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
                imp_occu = np.sum(denmat[0].diagonal()[0:v_norb]).real
                res.append((num_electron, eigval, imp_occu))
                if rank == 0:
                    print(num_electron, eigval, imp_occu, file=flog, flush=True)

            if nminus_direction:
                num_electron = nminus_list[i]
                if rank == 0:
                    write_fock_dec_by_N(ntot_v, num_electron, "fock_i.in")
                comm.Barrier()
                ed_fsolver(fcomm, rank, size)
                comm.Barrier()
                data = np.loadtxt('eigvals.dat', ndmin=2)
                eigval = data[0, 1]
                if eigval > eval_gs:
                    nminus_direction = False
                else:
                    nplus_direction = False
                    eval_gs = eigval
                    noccu_gs = num_electron
                data = np.loadtxt('denmat.dat', ndmin=2)
                tmp = (1, ntot_v, ntot_v)
                denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
                imp_occu = np.sum(denmat[0].diagonal()[0:v_norb]).real
                res.append((num_electron, eigval, imp_occu))
                if rank == 0:
                    print(num_electron, eigval, imp_occu, file=flog, flush=True)
        if rank == 0:
            flog.close()
            res.sort(key=lambda x: x[1])
            f = open('search_result.dat', 'w')
            for item in res:
                f.write("{:10d}{:20.10f}{:20.10f}\n".format(item[0], item[1], item[2]))
            f.close()
            print("edrixs >>> do_ed=0, Perform ED at occupancy: ", noccu_gs,
                  "with more accuracy", flush=True)
        # Do ED for the occupancy of ground state with more accuracy
        eval_shift = c_level * c_norb / noccu_gs
        emat_i[0:ntot_v, 0:ntot_v] += np.eye(ntot_v) * eval_shift
        emat_n[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level
        if rank == 0:
            write_emat(emat_i, 'hopping_i.in')
            write_emat(emat_n, 'hopping_n.in')
            write_config(
                ed_solver=ed_solver, num_val_orbs=ntot_v, neval=neval, nvector=nvector, ncv=ncv,
                idump=idump, maxiter=maxiter, min_ndim=min_ndim, eigval_tol=eigval_tol
            )
            write_fock_dec_by_N(ntot_v, noccu_gs, "fock_i.in")
        comm.Barrier()
        ed_fsolver(fcomm, rank, size)
        comm.Barrier()
        data = np.loadtxt('eigvals.dat', ndmin=2)
        eval_i = np.zeros(neval, dtype=float)
        eval_i[0:neval] = data[0:neval, 1]
        data = np.loadtxt('denmat.dat', ndmin=2)
        tmp = (nvector, ntot_v, ntot_v)
        denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)
        return eval_i, denmat, noccu_gs
    else:
        raise Exception("Unknown case of do_ed ", do_ed)


def xas_siam_fort(comm, shell_name, nbath, ominc, *, gamma_c=0.1,
                  v_noccu=1, thin=1.0, phi=0, pol_type=None,
                  num_gs=1, nkryl=200, temperature=1.0,
                  loc_axis=None, scatter_axis=None):
    """
    Calculate XAS for single impurity Anderson model (SIAM) with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') may indicate a :math:`L_3` edge transition from
        core :math:`2p_{3/2}` shell to valence :math:`3d` shell for Ni.
    nbath: int
        Number of bath sites.
    ominc: 1d float array
        Incident energy of photon.
    gamma_c: a float number or a 1d float array with the same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    v_noccu: int
        Total occupancy of valence shells.
    thin: float number
        The incident angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of the local scattering axis: scatter_axis[:,0].
    pol_type: list of tuples
        Type of polarization, options can be:

        - ('linear', alpha), linear polarization, where alpha is the angle between the
          polarization vector and the scattering plane in radians.

        - ('left', 0), left circular polarization.

        - ('right', 0), right circular polarization.

        - ('isotropic', 0). isotropic polarization.

        It will set pol_type=[('isotropic', 0)] if not provided.
    num_gs: int
        Number of initial states used in XAS calculations.
    nkryl: int
        Maximum number of poles obtained.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    xas: 2d array, shape=(len(ominc), len(pol_type))
        The calculated XAS spectra. The first dimension is for ominc, and the second dimension
        if for different polarizations.
    poles: list of dict, shape=(len(pol_type), )
        The calculated XAS poles for different polarizations.
    """
    _warn_legacy(
        'xas_siam_fort',
        'edrixs.models.model_siam + edrixs.solvers.get_ops/ed/xas',
    )
    from .fedrixs import xas_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']

    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    info_shell = info_atomic_shell()
    v_norb = info_shell[v_name][1]
    c_norb = info_shell[c_name][1]

    ntot_v = v_norb * (nbath + 1)
    ntot = ntot_v + c_norb
    if pol_type is None:
        pol_type = [('isotropic', 0)]
    if loc_axis is None:
        loc_axis = np.eye(3)
    else:
        loc_axis = np.array(loc_axis)
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    if rank == 0:
        print("edrixs >>> Running XAS ...", flush=True)
        write_config(num_val_orbs=ntot_v, num_core_orbs=c_norb,
                     num_gs=num_gs, nkryl=nkryl)
        write_fock_dec_by_N(ntot_v, v_noccu, "fock_i.in")
        write_fock_dec_by_N(ntot_v, v_noccu + 1, "fock_n.in")

    case = v_name + c_name
    tmp = get_trans_oper(case)
    npol, n, m = tmp.shape
    tmp_g = np.zeros((npol, n, m), dtype=complex)
    trans_mat = np.zeros((npol, ntot, ntot), dtype=complex)
    # Transform the transition operators to global-xyz axis
    # dipolar transition
    if npol == 3:
        for i in range(3):
            for j in range(3):
                tmp_g[i] += loc_axis[i, j] * tmp[j]
    # quadrupolar transition
    elif npol == 5:
        alpha, beta, gamma = rmat_to_euler(loc_axis)
        wignerD = get_wigner_dmat(4, alpha, beta, gamma)
        rotmat = np.dot(np.dot(tmat_r2c('d'), wignerD), np.conj(np.transpose(tmat_r2c('d'))))
        for i in range(5):
            for j in range(5):
                tmp_g[i] += rotmat[i, j] * tmp[j]
    else:
        raise Exception("Have NOT implemented this case: ", npol)
    trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

    n_om = len(ominc)
    gamma_core = np.zeros(n_om, dtype=float)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_om) * gamma_c
    else:
        gamma_core[:] = gamma_c

    # loop over different polarization
    xas = np.zeros((n_om, len(pol_type)), dtype=float)
    poles = []
    comm.Barrier()
    for it, (pt, alpha) in enumerate(pol_type):
        if pt.strip() == 'left' or pt.strip() == 'right' or pt.strip() == 'linear':
            if rank == 0:
                print("edrixs >>> Loop over for polarization: ", it, pt, flush=True)
                kvec = unit_wavevector(thin, phi, scatter_axis, 'in')
                polvec = np.zeros(npol, dtype=complex)
                pol = dipole_polvec_xas(thin, phi, alpha, scatter_axis, pt)
                if npol == 3:  # Dipolar transition
                    polvec[:] = pol
                if npol == 5:  # Quadrupolar transition
                    polvec[:] = quadrupole_polvec(pol, kvec)

                trans = np.zeros((ntot, ntot), dtype=complex)
                for i in range(npol):
                    trans[:, :] += trans_mat[i] * polvec[i]
                write_emat(trans, 'transop_xas.in')

            # call XAS solver in fedrixs
            comm.Barrier()
            xas_fsolver(fcomm, rank, size)
            comm.Barrier()

            file_list = ['xas_poles.' + str(i+1) for i in range(num_gs)]
            pole_dict = read_poles_from_file(file_list)
            poles.append(pole_dict)
            xas[:, it] = get_spectra_from_poles(pole_dict, ominc, gamma_core, temperature)
        elif pt.strip() == 'isotropic':
            pole_dicts = []
            for k in range(npol):
                if rank == 0:
                    print("edrixs >>> Loop over for polarization: ", it, pt, flush=True)
                    print("edrixs >>> Isotropic, component: ", k, flush=True)
                    write_emat(trans_mat[k], 'transop_xas.in')

                # call XAS solver in fedrixs
                comm.Barrier()
                xas_fsolver(fcomm, rank, size)
                comm.Barrier()

                file_list = ['xas_poles.' + str(i+1) for i in range(num_gs)]
                pole_tmp = read_poles_from_file(file_list)
                xas[:, it] += get_spectra_from_poles(pole_tmp, ominc, gamma_core, temperature)
                pole_dicts.append(pole_tmp)
            xas[:, it] = xas[:, it] / npol
            poles.append(merge_pole_dicts(pole_dicts))
        else:
            raise Exception("Unknown polarization type: ", pt)

    return xas, poles


def rixs_siam_fort(comm, shell_name, nbath, ominc, eloss, *, gamma_c=0.1, gamma_f=0.1,
                   v_noccu=1, thin=1.0, thout=1.0, phi=0, pol_type=None, num_gs=1,
                   nkryl=200, linsys_max=1000, linsys_tol=1e-10, temperature=1.0,
                   loc_axis=None, scatter_axis=None):
    """
    Calculate RIXS for single impurity Anderson model with Fortran solver.

    Parameters
    ----------
    comm: MPI_comm
        MPI communicator.
    shell_name: tuple of two strings
        Names of valence and core shells. The 1st (2nd) string in the tuple is for the
        valence (core) shell.

        - The 1st string can only be 's', 'p', 't2g', 'd', 'f',

        - The 2nd string can be 's', 'p', 'p12', 'p32', 'd', 'd32', 'd52',
          'f', 'f52', 'f72'.

        For example: shell_name=('d', 'p32') may indicate a :math:`L_3` edge transition from
        core :math:`2p_{3/2}` shell to valence :math:`3d` shell for Ni.
    nbath: int
        Number of bath sites.
    ominc: 1d float array
        Incident energy of photon.
    eloss: 1d float array
        Energy loss.
    gamma_c: a float number or a 1d float array with same shape as ominc.
        The core-hole life-time broadening factor. It can be a constant value
        or incident energy dependent.
    gamma_f: a float number or a 1d float array with same shape as eloss.
        The final states life-time broadening factor. It can be a constant value
        or energy loss dependent.
    v_noccu: int
        Total occupancy of valence shells.
    thin: float number
        The incident angle of photon (in radian).
    thout: float number
        The scattered angle of photon (in radian).
    phi: float number
        Azimuthal angle (in radian), defined with respect to the
        :math:`x`-axis of scattering axis: scatter_axis[:,0].
    pol_type: list of 4-elements-tuples
        Type of polarizations. It has the following form:

        (str1, alpha, str2, beta)

        where, str1 (str2) can be 'linear', 'left', 'right', and alpha (beta) is
        the angle (in radians) between the linear polarization vector and the scattering plane.

        It will set pol_type=[('linear', 0, 'linear', 0)] if not provided.
    num_gs: int
        Number of initial states used in RIXS calculations.
    nkryl: int
        Maximum number of poles obtained.
    linsys_max: int
        Maximum iterations of solving linear equations.
    linsys_tol: float
        Convergence for solving linear equations.
    temperature: float number
        Temperature (in K) for boltzmann distribution.
    loc_axis: 3*3 float array
        The local axis with respect to which local orbitals are defined.

        - x: local_axis[:,0],

        - y: local_axis[:,1],

        - z: local_axis[:,2].

        It will be an identity matrix if not provided.
    scatter_axis: 3*3 float array
        The local axis defining the scattering geometry. The scattering plane is defined in
        the local :math:`zx`-plane.

        - local :math:`x`-axis: scatter_axis[:,0]

        - local :math:`y`-axis: scatter_axis[:,1]

        - local :math:`z`-axis: scatter_axis[:,2]

        It will be set to an identity matrix if not provided.

    Returns
    -------
    rixs: 3d float array, shape=(len(ominc), len(eloss), len(pol_type))
        The calculated RIXS spectra. The 1st dimension is for the incident energy,
        the 2nd dimension is for the energy loss and the 3rd dimension is for
        different polarizations.
    poles: 2d list of dict, shape=(len(ominc), len(pol_type))
        The calculated RIXS poles. The 1st dimension is for incident energy, and the
        2nd dimension is for different polarizations.
    """
    _warn_legacy(
        'rixs_siam_fort',
        'edrixs.models.model_siam + edrixs.solvers.get_ops/ed/rixs',
    )
    from .fedrixs import rixs_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']
    v_name = shell_name[0].strip()
    c_name = shell_name[1].strip()
    if v_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    info_shell = info_atomic_shell()
    v_norb = info_shell[v_name][1]
    c_norb = info_shell[c_name][1]
    ntot_v = v_norb * (nbath + 1)
    ntot = ntot_v + c_norb

    if pol_type is None:
        pol_type = [('linear', 0, 'linear', 0)]
    if loc_axis is None:
        loc_axis = np.eye(3)
    else:
        loc_axis = np.array(loc_axis)
    if scatter_axis is None:
        scatter_axis = np.eye(3)
    else:
        scatter_axis = np.array(scatter_axis)

    if rank == 0:
        print("edrixs >>> Running RIXS ...", flush=True)
        write_fock_dec_by_N(ntot_v, v_noccu, "fock_i.in")
        write_fock_dec_by_N(ntot_v, v_noccu + 1, "fock_n.in")
        write_fock_dec_by_N(ntot_v, v_noccu, "fock_f.in")

        case = v_name + c_name
        tmp = get_trans_oper(case)
        npol, n, m = tmp.shape
        tmp_g = np.zeros((npol, n, m), dtype=complex)
        trans_mat = np.zeros((npol, ntot, ntot), dtype=complex)
        # Transform the transition operators to global-xyz axis
        # dipolar transition
        if npol == 3:
            for i in range(3):
                for j in range(3):
                    tmp_g[i] += loc_axis[i, j] * tmp[j]
        # quadrupolar transition
        elif npol == 5:
            alpha, beta, gamma = rmat_to_euler(loc_axis)
            wignerD = get_wigner_dmat(4, alpha, beta, gamma)
            rotmat = np.dot(np.dot(tmat_r2c('d'), wignerD), np.conj(np.transpose(tmat_r2c('d'))))
            for i in range(5):
                for j in range(5):
                    tmp_g[i] += rotmat[i, j] * tmp[j]
        else:
            raise Exception("Have NOT implemented this case: ", npol)
        trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

    n_om = len(ominc)
    neloss = len(eloss)
    gamma_core = np.zeros(n_om, dtype=float)
    if np.isscalar(gamma_c):
        gamma_core[:] = np.ones(n_om) * gamma_c
    else:
        gamma_core[:] = gamma_c
    gamma_final = np.zeros(neloss, dtype=float)
    if np.isscalar(gamma_f):
        gamma_final[:] = np.ones(neloss) * gamma_f
    else:
        gamma_final[:] = gamma_f

    # loop over different polarization
    rixs = np.zeros((n_om, neloss, len(pol_type)), dtype=float)
    poles = []
    comm.Barrier()
    # loop over different polarization
    for iom, omega in enumerate(ominc):
        if rank == 0:
            write_config(
                num_val_orbs=ntot_v, num_core_orbs=c_norb,
                omega_in=omega, gamma_in=gamma_core[iom],
                num_gs=num_gs, nkryl=nkryl, linsys_max=linsys_max,
                linsys_tol=linsys_tol
            )
        poles_per_om = []
        # loop over polarization
        for ip, (it, alpha, jt, beta) in enumerate(pol_type):
            if rank == 0:
                print(flush=True)
                print("edrixs >>> Calculate RIXS for incident energy: ", omega, flush=True)
                print("edrixs >>> Polarization: ", ip, flush=True)
                polvec_i = np.zeros(npol, dtype=complex)
                polvec_f = np.zeros(npol, dtype=complex)
                ei, ef = dipole_polvec_rixs(thin, thout, phi, alpha, beta,
                                            scatter_axis, (it, jt))
                # dipolar transition
                if npol == 3:
                    polvec_i[:] = ei
                    polvec_f[:] = ef
                # quadrupolar transition
                elif npol == 5:
                    ki = unit_wavevector(thin, phi, scatter_axis, direction='in')
                    kf = unit_wavevector(thout, phi, scatter_axis, direction='out')
                    polvec_i[:] = quadrupole_polvec(ei, ki)
                    polvec_f[:] = quadrupole_polvec(ef, kf)
                else:
                    raise Exception("Have NOT implemented this type of transition operators")
                trans_i = np.zeros((ntot, ntot), dtype=complex)
                trans_f = np.zeros((ntot, ntot), dtype=complex)
                for i in range(npol):
                    trans_i[:, :] += trans_mat[i] * polvec_i[i]
                write_emat(trans_i, 'transop_rixs_i.in')
                for i in range(npol):
                    trans_f[:, :] += trans_mat[i] * polvec_f[i]
                write_emat(np.conj(np.transpose(trans_f)), 'transop_rixs_f.in')

            # call RIXS solver in fedrixs
            comm.Barrier()
            rixs_fsolver(fcomm, rank, size)
            comm.Barrier()

            file_list = ['rixs_poles.' + str(i+1) for i in range(num_gs)]
            pole_dict = read_poles_from_file(file_list)
            poles_per_om.append(pole_dict)
            rixs[iom, :, ip] = get_spectra_from_poles(pole_dict, eloss,
                                                      gamma_final, temperature)

        poles.append(poles_per_om)

    return rixs, poles
