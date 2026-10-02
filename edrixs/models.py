"""Backend-independent physical model construction for EDRIXS.

The model functions in this module describe orbital-space Hamiltonians,
interactions, Fock-basis specifications, and photon-transition matrices without
constructing backend-owned many-body operators.
"""

from __future__ import annotations

import numpy as np

from .angular_momentum import get_lx, get_ly, get_lz, get_sx, get_sy, get_sz
from .angular_momentum import get_wigner_dmat, rmat_to_euler
from .basis_transform import cb_op, tmat_c2r, tmat_r2c
from .coulomb_utensor import get_F0, get_umat_slater, get_umat_slater_3shells
from .fock_basis import FockBasisSpec
from .photon_transition import get_trans_oper
from .soc import atom_hsoc
from .utils import (CT_imp_bath, CT_imp_bath_core_hole, info_atomic_shell,
                    slater_integrals_name)
from ._solvers_helpers import (
    _embed_impurity_core_umat,
    _embed_impurity_core_umat_sparse,
    _rotated_transition_blocks,
    _umat_dense_to_sparse,
    _valence_zeeman_matrix,
)

__all__ = ['model_1v1c', 'model_2v1c', 'model_siam', 'model_siam_2d1p']


def _print_slater_summary(slater_name, slater_i, slater_n):
    """Print the initial/intermediate Slater-integral table."""
    print()
    print("    Summary of Slater integrals:")
    print("    ------------------------------")
    print("    Terms,   Initial Hamiltonian,  Intermediate Hamiltonian")
    for name, s_i, s_n in zip(slater_name, slater_i, slater_n):
        print("    ", name, ":  {:20.10f}{:20.10f}".format(s_i, s_n))
    print()


def model_1v1c(shell_name, *, shell_level=None, v_soc=None, c_soc=0,
               v_noccu=1, slater=None, ext_B=None, on_which='spin',
               v_cfmat=None, v_othermat=None, loc_axis=None, verbose=False,
               sparse_U=False, tol=1E-10):
    """
    Set up orbital-space data and Fock-basis metadata for a 1v1c problem.

    This routine defines the physical one-valence-shell/one-core-shell problem
    independently of the numerical backend. It constructs the one-body orbital
    matrices, Coulomb tensors, Fock-basis specifications, and orbital-space
    transition matrices, but it does not build many-body Hamiltonians and does not
    diagonalize anything.

    Parameters
    ----------
    shell_name : tuple of str
        Names of the valence and core shells.

    shell_level, v_soc, c_soc, v_noccu, slater, ext_B, on_which,
    v_cfmat, v_othermat, loc_axis
        Same physical model parameters as ed_1v1c_py.

    verbose : bool, optional
        If True, print a setup summary (Slater integrals and Hilbert-space
        dimensions). Default False.

    sparse_U : bool, optional
        If False, return dense rank-4 Coulomb tensors. If True, return each
        Coulomb tensor as a sparse matrix with shape (ntot*ntot, ntot*ntot),
        using the flattened convention
        row = lorb * ntot + korb and col = jorb * ntot + iorb.

    tol : float, optional
        Threshold used when converting dense Coulomb tensors to sparse format.

    Returns
    -------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, shift
        Backend-independent problem definition. ``basis_i`` and ``basis_n`` are
        compact Fock-basis specifications with a filled initial core and one
        intermediate core hole. Occupancy arguments count valence electrons
        only. Both one-body matrices have shape (ntot, ntot), both dense
        Coulomb tensors have shape (ntot,) * 4, and trans_mat has shape
        (npol, ntot, ntot). Initial core interactions are included.
        ``shift`` is the filled core shell-level energy,
        ``c_norb * shell_level[1]`` (zero if unspecified), excluding Coulomb
        contributions.
        Pass the first seven outputs to :func:`~edrixs.solvers.get_ops` and
        ``shift`` to dense, SciPy, or PETSc :func:`~edrixs.solvers.ed`.
        The returned integrals retain the full, unshifted energies.
    """
    if verbose:
        print("edrixs >>> Setting up 1v1c problem ...")

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52',
                      'f', 'f52', 'f72']

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
    ntot = v_norb + c_norb

    # Both sectors explicitly include the core orbitals.
    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)

    # Coulomb interaction.
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

    if verbose:
        _print_slater_summary(slater_name, slater_i, slater_n)

    case = v_name + c_name
    umat_i = get_umat_slater(case, *slater_i)
    umat_n = get_umat_slater(case, *slater_n)

    if sparse_U:
        umat_i = _umat_dense_to_sparse(umat_i, tol=tol)
        umat_n = _umat_dense_to_sparse(umat_n, tol=tol)

    # Spin-orbit coupling.
    if v_soc is not None:
        emat_i[0:v_norb, 0:v_norb] += atom_hsoc(v_name, v_soc[0])
        emat_n[0:v_norb, 0:v_norb] += atom_hsoc(v_name, v_soc[1])

    # For split core shells such as p12/p32, d32/d52, f52/f72, the SOC is
    # already encoded by the shell choice.
    if c_name in ['p', 'd', 'f']:
        emat_i[v_norb:ntot, v_norb:ntot] += atom_hsoc(c_name, c_soc)
        emat_n[v_norb:ntot, v_norb:ntot] += atom_hsoc(c_name, c_soc)

    # Crystal field and additional one-body terms.
    if v_cfmat is not None:
        emat_i[0:v_norb, 0:v_norb] += np.asarray(v_cfmat)
        emat_n[0:v_norb, 0:v_norb] += np.asarray(v_cfmat)

    if v_othermat is not None:
        emat_i[0:v_norb, 0:v_norb] += np.asarray(v_othermat)
        emat_n[0:v_norb, 0:v_norb] += np.asarray(v_othermat)

    # Shell levels.
    if shell_level is not None:
        emat_i[0:v_norb, 0:v_norb] += np.eye(v_norb) * shell_level[0]

        emat_n[0:v_norb, 0:v_norb] += np.eye(v_norb) * shell_level[0]
        emat_i[v_norb:ntot, v_norb:ntot] += np.eye(c_norb) * shell_level[1]
        emat_n[v_norb:ntot, v_norb:ntot] += np.eye(c_norb) * shell_level[1]

    # External field on valence shell.
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
            zeeman = (
                ext_B[0] * (lx + 2 * sx)
                + ext_B[1] * (ly + 2 * sy)
                + ext_B[2] * (lz + 2 * sz)
            )
        else:
            raise Exception("Unknown value of on_which", on_which)

        emat_i[0:v_norb, 0:v_norb] += zeeman
        emat_n[0:v_norb, 0:v_norb] += zeeman

    # Fock-basis metadata.
    basis_i = FockBasisSpec.from_args(v_norb, v_noccu, c_norb, c_norb)
    basis_n = FockBasisSpec.from_args(v_norb, v_noccu + 1, c_norb, c_norb - 1)

    if verbose:
        print("edrixs >>> Dimension of the initial Hamiltonian: ", len(basis_i))
        print("edrixs >>> Dimension of the intermediate Hamiltonian: ", len(basis_n))

    # Transition operators in local coordinates, rotated to global coordinates.
    if loc_axis is not None:
        local_axis = np.asarray(loc_axis)
    else:
        local_axis = np.eye(3)

    tmp = get_trans_oper(case)
    npol, n, m = tmp.shape
    tmp_g = np.zeros((npol, n, m), dtype=complex)

    if npol == 3:
        for i in range(3):
            for j in range(3):
                tmp_g[i] += local_axis[i, j] * tmp[j]
    elif npol == 5:
        alpha, beta, gamma = rmat_to_euler(local_axis)
        wignerD = get_wigner_dmat(4, alpha, beta, gamma)
        rotmat = np.dot(
            np.dot(tmat_r2c('d'), wignerD),
            np.conj(np.transpose(tmat_r2c('d')))
        )
        for i in range(5):
            for j in range(5):
                tmp_g[i] += rotmat[i, j] * tmp[j]
    else:
        raise Exception("Have NOT implemented this case: ", npol)

    trans_mat = np.zeros((npol, ntot, ntot), dtype=complex)
    for i in range(npol):
        trans_mat[i, 0:v_norb, v_norb:ntot] = tmp_g[i]

    if verbose:
        print("edrixs >>> 1v1c setup Done !")

    shift = 0.0 if shell_level is None else float(c_norb * shell_level[1])
    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, shift


def model_2v1c(
    shell_name, *, shell_level=None,
    v1_soc=None, v2_soc=None, c_soc=0, v_tot_noccu=1, slater=None,
    v1_ext_B=None, v2_ext_B=None, v1_on_which='spin',
    v2_on_which='spin', v1_cfmat=None, v2_cfmat=None,
    v1_othermat=None, v2_othermat=None, hopping_v1v2=None,
    trans_to_which=1, loc_axis=None, verbose=False, sparse_U=False, tol=1E-10
):
    """
    Set up orbital-space data and Fock-basis metadata for a 2v1c problem.

    This is the backend-neutral setup analogue of the 2-valence-shell,
    1-core-shell Fortran ED/RIXS input construction.  It does not build
    many-body Hamiltonians and does not diagonalize anything. If sparse_U=True,
    the returned Coulomb tensors are sparse flattened matrices rather than
    dense rank-4 arrays.

    Parameters
    ----------
    shell_name : tuple of str
        Names of the first valence shell, second valence shell, and core shell,
        in that order.  Each valence shell must be one of ``'s'``, ``'p'``,
        ``'t2g'``, ``'d'``, or ``'f'``.  The core shell may additionally be a
        spin-orbit-split shell such as ``'p12'``, ``'p32'``, ``'d32'``,
        ``'d52'``, ``'f52'``, or ``'f72'``.
    shell_level : sequence of three float or None, optional
        One-body energy levels of the first valence, second valence, and core
        shells.  All levels default to zero.
    v1_soc, v2_soc : sequence of two float or None, optional
        Spin-orbit coupling strengths for the corresponding valence shell in
        the initial and intermediate Hamiltonians, respectively.  Spin-orbit
        coupling is omitted when the value is ``None``.
    c_soc : float, optional
        Core-shell spin-orbit coupling strength.  This is applied to unsplit
        ``'p'``, ``'d'``, and ``'f'`` core shells.  Default is zero.
    v_tot_noccu : int, optional
        Total number of electrons occupying the two valence shells in the
        initial state.  Default is one.
    slater : tuple of two sequences of float or None, optional
        Slater integrals for the initial and intermediate Hamiltonians.  Each
        sequence follows the order returned by
        ``slater_integrals_name(shell_name, ('v1', 'v2', 'c1'))``.  Missing
        trailing values and all values when ``slater`` is ``None`` are treated
        as zero.
    v1_ext_B, v2_ext_B : array-like of float, shape (3,), or None, optional
        Zeeman/exchange-field coefficients in the global Cartesian frame for
        the corresponding valence shell.  No field is applied when the value
        is ``None``.
    v1_on_which, v2_on_which : {'spin', 'orbital', 'both'}, optional
        Sector on which the corresponding external field acts.  Default is
        ``'spin'``.
    v1_cfmat, v2_cfmat : array-like of complex or None, optional
        Crystal-field matrix for the corresponding valence shell.  Its shape
        must match that shell's number of spin-orbitals.  Defaults to zero.
    v1_othermat, v2_othermat : array-like of complex or None, optional
        Additional one-body matrix for the corresponding valence shell.  Its
        shape must match that shell's number of spin-orbitals.  Defaults to
        zero.
    hopping_v1v2 : array-like of complex or None, optional
        Hopping matrix from the second valence shell to the first, with shape
        ``(v1_norb, v2_norb)``.  Its Hermitian conjugate is inserted
        automatically.  Defaults to zero.
    trans_to_which : {1, 2}, optional
        Valence shell reached by the core-to-valence transition.  Default is
        the first valence shell.
    loc_axis : array-like of float, shape (3, 3), or None, optional
        Rotation from local crystal-field coordinates to the global frame for
        the transition operators.  The identity rotation is used by default.
    verbose : bool, optional
        If True, print a setup summary (Slater integrals and Hilbert-space
        dimensions).  Default is False.
    sparse_U : bool, optional
        If False, return dense rank-4 Coulomb tensors.  If True, return CSR
        matrices using the flattened convention
        ``row = lorb * norb + korb`` and ``col = jorb * norb + iorb``.
        Default is False.
    tol : float, optional
        Absolute-value threshold used when converting Coulomb tensors to
        sparse form.  Default is ``1e-10``.

    Returns
    -------
    emat_i : numpy.ndarray
        Initial-state one-body matrix over both valence shells and the core,
        with shape ``(ntot, ntot)``.
    umat_i : numpy.ndarray or scipy.sparse.csr_matrix
        Initial-state Coulomb interaction including core interactions. The dense
        representation has shape ``(ntot,) * 4``; the sparse
        representation has shape ``(ntot**2, ntot**2)``.
    basis_i : FockBasisSpec
        Initial-state Fock-basis specification with ``v_tot_noccu`` valence
        electrons and a fully occupied core shell.
    emat_n : numpy.ndarray
        Intermediate-state one-body matrix over both valence shells and the
        core shell, with shape ``(ntot, ntot)``.
    umat_n : numpy.ndarray or scipy.sparse.csr_matrix
        Intermediate-state Coulomb interaction.  The dense representation has
        shape ``(ntot,) * 4``; the sparse representation has shape
        ``(ntot**2, ntot**2)``.
    basis_n : FockBasisSpec
        Intermediate-state Fock-basis specification with one additional
        valence electron and one core hole.
    trans_mat : numpy.ndarray
        Core-to-valence transition matrices in the global frame, with shape
        ``(npol, ntot, ntot)``.

    shift : float
        Filled core shell-level energy, ``c_norb * shell_level[2]`` (zero if
        unspecified), excluding Coulomb contributions. Pass the first seven
        outputs to :func:`~edrixs.solvers.get_ops` and ``shift`` to dense, SciPy, or
        PETSc :func:`~edrixs.solvers.ed`. Integrals remain unshifted.

    See Also
    --------
    edrixs.models.model_1v1c
        Construct a model with one valence shell and one core shell.
    edrixs.models.model_siam
        Construct a single-impurity Anderson model.
    edrixs.solvers.get_ops
        Convert a model definition into backend-owned many-body operators.

    Examples
    --------
    Construct a model for two valence shells and build its SciPy operators:

    >>> import edrixs
    >>> *problem, shift = edrixs.model_2v1c(
    ...     ('d', 'p', 's'), v_tot_noccu=2, trans_to_which=2
    ... )
    >>> hmat_i, hmat_n, trans_ops = edrixs.get_ops(
    ...     *problem, backend='scipy'
    ... )
    """
    if verbose:
        print("edrixs >>> Setting up 2v1c problem ...")

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = [
        's', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52',
        'f', 'f52', 'f72'
    ]

    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()

    if v1_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v1_name)
    if v2_name not in v_name_options:
        raise Exception("NOT supported type of valence shell: ", v2_name)
    if c_name not in c_name_options:
        raise Exception("NOT supported type of core shell: ", c_name)

    info_shell = info_atomic_shell()

    v1_orbl = info_shell[v1_name][0]
    v2_orbl = info_shell[v2_name][0]

    v1_norb = info_shell[v1_name][1]
    v2_norb = info_shell[v2_name][1]
    c_norb = info_shell[c_name][1]

    v1v2_norb = v1_norb + v2_norb
    ntot = v1v2_norb + c_norb

    slater_name = slater_integrals_name(
        (v1_name, v2_name, c_name), ('v1', 'v2', 'c1')
    )
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

    if verbose:
        _print_slater_summary(slater_name, slater_i, slater_n)

    umat_i = get_umat_slater_3shells(
        (v1_name, v2_name, c_name), *slater_i
    )
    umat_n = get_umat_slater_3shells(
        (v1_name, v2_name, c_name), *slater_n
    )

    if sparse_U:
        umat_i = _umat_dense_to_sparse(umat_i, tol=tol)
        umat_n = _umat_dense_to_sparse(umat_n, tol=tol)

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)

    # Spin-orbit coupling.
    if v1_soc is not None and v1_name in ['p', 'd', 't2g', 'f']:
        emat_i[0:v1_norb, 0:v1_norb] += atom_hsoc(v1_name, v1_soc[0])
        emat_n[0:v1_norb, 0:v1_norb] += atom_hsoc(v1_name, v1_soc[1])

    if v2_soc is not None and v2_name in ['p', 'd', 't2g', 'f']:
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += atom_hsoc(
            v2_name, v2_soc[0]
        )
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += atom_hsoc(
            v2_name, v2_soc[1]
        )

    if c_name in ['p', 'd', 'f']:
        emat_i[v1v2_norb:ntot, v1v2_norb:ntot] += atom_hsoc(
            c_name, c_soc
        )
        emat_n[v1v2_norb:ntot, v1v2_norb:ntot] += atom_hsoc(
            c_name, c_soc
        )

    # Crystal fields.
    if v1_cfmat is not None:
        emat_i[0:v1_norb, 0:v1_norb] += np.asarray(v1_cfmat)
        emat_n[0:v1_norb, 0:v1_norb] += np.asarray(v1_cfmat)

    if v2_cfmat is not None:
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.asarray(
            v2_cfmat
        )
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.asarray(
            v2_cfmat
        )

    # Other one-body terms.
    if v1_othermat is not None:
        emat_i[0:v1_norb, 0:v1_norb] += np.asarray(v1_othermat)
        emat_n[0:v1_norb, 0:v1_norb] += np.asarray(v1_othermat)

    if v2_othermat is not None:
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.asarray(
            v2_othermat
        )
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.asarray(
            v2_othermat
        )

    # Shell levels, including the explicit filled core.
    if shell_level is not None:
        emat_i[0:v1_norb, 0:v1_norb] += np.eye(v1_norb) * shell_level[0]
        emat_n[0:v1_norb, 0:v1_norb] += np.eye(v1_norb) * shell_level[0]

        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += (
            np.eye(v2_norb) * shell_level[1]
        )
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += (
            np.eye(v2_norb) * shell_level[1]
        )

        emat_i[v1v2_norb:ntot, v1v2_norb:ntot] += (
            np.eye(c_norb) * shell_level[2]
        )
        emat_n[v1v2_norb:ntot, v1v2_norb:ntot] += (
            np.eye(c_norb) * shell_level[2]
        )

    # Zeeman fields.
    for name, lval, ext_B, which, i1, i2 in [
        (v1_name, v1_orbl, v1_ext_B, v1_on_which, 0, v1_norb),
        (
            v2_name, v2_orbl, v2_ext_B, v2_on_which,
            v1_norb, v1v2_norb
        )
    ]:
        if ext_B is None:
            continue

        zeeman = _valence_zeeman_matrix(name, lval, ext_B, which)
        emat_i[i1:i2, i1:i2] += zeeman
        emat_n[i1:i2, i1:i2] += zeeman

    # Hopping between the two valence shells.
    if hopping_v1v2 is not None:
        hopping_v1v2 = np.asarray(hopping_v1v2)
        emat_i[0:v1_norb, v1_norb:v1v2_norb] += hopping_v1v2
        emat_i[v1_norb:v1v2_norb, 0:v1_norb] += np.conj(
            np.transpose(hopping_v1v2)
        )

        emat_n[0:v1_norb, v1_norb:v1v2_norb] += hopping_v1v2
        emat_n[v1_norb:v1v2_norb, 0:v1_norb] += np.conj(
            np.transpose(hopping_v1v2)
        )

    basis_i = FockBasisSpec.from_args(v1v2_norb, v_tot_noccu, c_norb, c_norb)
    basis_n = FockBasisSpec.from_args(
        v1v2_norb, v_tot_noccu + 1, c_norb, c_norb - 1
    )

    if verbose:
        print("edrixs >>> Dimension of the initial Hamiltonian: ", len(basis_i))
        print("edrixs >>> Dimension of the intermediate Hamiltonian: ", len(basis_n))

    if trans_to_which == 1:
        case = v1_name + c_name
        tmp_g = _rotated_transition_blocks(case, loc_axis)
        trans_mat = np.zeros((tmp_g.shape[0], ntot, ntot), dtype=complex)
        trans_mat[:, 0:v1_norb, v1v2_norb:ntot] = tmp_g
    elif trans_to_which == 2:
        case = v2_name + c_name
        tmp_g = _rotated_transition_blocks(case, loc_axis)
        trans_mat = np.zeros((tmp_g.shape[0], ntot, ntot), dtype=complex)
        trans_mat[:, v1_norb:v1v2_norb, v1v2_norb:ntot] = tmp_g
    else:
        raise Exception("Unknown trans_to_which: ", trans_to_which)

    if verbose:
        print("edrixs >>> 2v1c setup Done !")

    shift = 0.0 if shell_level is None else float(c_norb * shell_level[2])
    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, shift


def model_siam(
    shell_name, nbath, *, siam_type=0, v_noccu=1, static_core_pot=0,
    c_level=0, c_soc=0, trans_c2n=None, imp_mat=None, imp_mat_n=None,
    bath_level=None, bath_level_n=None, hyb=None, hyb_n=None,
    hopping=None, hopping_n=None, slater=None, ext_B=None,
    on_which='spin', loc_axis=None, verbose=False, sparse_U=False, tol=1E-10
):
    """
    Set up a single-impurity Anderson model (SIAM) for ED, XAS, and RIXS.

    Construct backend-independent one-body matrices, Coulomb interactions,
    Fock-basis specifications, and core-to-impurity transition matrices. This
    is the setup analogue of :func:`~edrixs.solvers.ed_siam_fort` at fixed
    ``v_noccu``; it does not search over occupancies, build many-body
    Hamiltonians, or diagonalize them. Pass the returned tuple to
    :func:`~edrixs.solvers.get_ops` to construct the many-body operators
    with the SciPy, dense, or PETSc backend.

    All energies, Slater integrals, hybridizations, spin-orbit coupling
    constants, and magnetic-field components are in eV. Let ``v_norb`` and
    ``c_norb`` be the numbers of impurity and core spin-orbitals, respectively,
    ``ntot_v = (nbath + 1) * v_norb``, and ``ntot = ntot_v + c_norb``.
    Orbitals are ordered as impurity, successive bath sites, then core.
    Each bath site has ``v_norb`` spin-orbitals and no Coulomb interaction.

    Parameters
    ----------
    shell_name : tuple of two str
        Impurity valence and core shell names, in that order. The valence
        shell can be 's', 'p', 't2g', 'd', or 'f'. The core shell can be
        's', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', or
        'f72'. For example, ``('d', 'p32')`` describes an L3-edge transition
        from a p3/2 core shell to a d impurity.
    nbath : int
        Number of bath sites.
    siam_type : {0, 1}, optional
        Parameterization of the impurity and bath one-body terms. For 0
        (default), use ``imp_mat``, ``bath_level``, and ``hyb`` to describe
        diagonal hybridization. For 1, use the full ``hopping`` matrix,
        allowing off-diagonal hybridization. The corresponding ``*_n``
        parameters specify intermediate-state terms. Parameters belonging
        to the other parameterization are ignored.
    v_noccu : int, optional
        Total number of electrons in the impurity and bath orbitals in the
        initial state, excluding the filled core. Default is 1. Must satisfy
        ``0 <= v_noccu < ntot_v`` so that the intermediate state can
        accommodate one additional electron.
    static_core_pot : float, optional
        Static core-hole potential, subtracted from every impurity diagonal
        element in the intermediate state. Positive values are attractive.
        This is added independently of the core-valence Slater interactions.
        Default is 0.
    c_level : float, optional
        Core-shell energy level, added to each core diagonal element in both
        states. Default is 0.
    c_soc : float, optional
        Core spin-orbit coupling strength, applied for unsplit 'p', 'd', and
        'f' core shells only. For split shells such as 'p32', the shell
        choice already selects the total-angular-momentum subspace.
        Default is 0.
    trans_c2n : array_like of complex, shape (v_norb, v_norb), optional
        Unitary transformation from complex spherical harmonics to the
        basis used for the impurity and bath one-body input parameters.
        Its columns are the input basis vectors expressed in spherical
        harmonics. The same transformation is used on the impurity and
        every bath site to convert the one-body matrices back to spherical
        harmonics. Default is the identity.
    imp_mat : array_like of complex, shape (v_norb, v_norb), optional
        Initial-state impurity one-body matrix for ``siam_type=0``, including
        the impurity energy level, crystal field, and valence spin-orbit
        coupling as needed. Default is zero.
    imp_mat_n : array_like of complex, shape (v_norb, v_norb), optional
        Intermediate-state impurity matrix for ``siam_type=0``. If None,
        use ``imp_mat``.
    bath_level : numpy.ndarray, shape (nbath, v_norb), optional
        Initial-state bath energy levels for ``siam_type=0``. The first
        index selects the bath site and the second its spin-orbital.
        Default is zero.
    bath_level_n : numpy.ndarray, shape (nbath, v_norb), optional
        Intermediate-state bath levels for ``siam_type=0``. If None, use
        ``bath_level``.
    hyb : numpy.ndarray, shape (nbath, v_norb), optional
        Initial-state hybridization amplitudes for ``siam_type=0``.
        ``hyb[i, j]`` is the matrix element from spin-orbital ``j`` of bath
        site ``i`` to impurity spin-orbital ``j``; its complex conjugate is
        inserted for the reverse hopping. Default is zero.
    hyb_n : numpy.ndarray, shape (nbath, v_norb), optional
        Intermediate-state hybridization amplitudes for ``siam_type=0``.
        If None, use ``hyb``.
    hopping : array_like of complex, shape (ntot_v, ntot_v), optional
        Full initial-state one-body matrix for ``siam_type=1``, including
        impurity terms, bath levels, and hybridization. Supply a Hermitian
        matrix in impurity-then-bath order. Default is zero.
    hopping_n : array_like of complex, shape (ntot_v, ntot_v), optional
        Full intermediate-state impurity and bath matrix for ``siam_type=1``.
        If None, use ``hopping``.
    slater : tuple of two lists, optional
        Slater integrals for the initial and intermediate states, respectively.
        Each list is ordered as ``[FX_vv, FX_vc, GX_vc, FX_cc]``, with ranks
        X in ascending order within each group; v denotes the impurity and
        c the core. Use ``slater_integrals_name(shell_name, ('v', 'c'))``
        to obtain the names in the required order. Trailing zeros may be
        omitted; entries beyond the required number are ignored. For
        example, d-impurity/p-core integrals begin with
        ``[F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp]``.
        Initial-state core-valence and core-core terms act on the explicit
        filled core. Default is all zeros.
    ext_B : array_like of float, shape (3,), optional
        External Zeeman/exchange field components along the global x, y,
        and z axes, applied only to the impurity in both states. Components
        are energy coefficients, not magnetic fields in tesla. Default is
        zero.
    on_which : {'spin', 'orbital', 'both'}, optional
        Couple ``ext_B`` to ``2S``, ``L``, or ``L + 2S``, respectively,
        with a positive sign in the Hamiltonian. Default is 'spin'.
    loc_axis : array_like of float, shape (3, 3), optional
        Rotation from local crystal-field coordinates to the global frame
        for the transition operators. Columns are the local axes expressed
        in global coordinates. Default is the identity rotation.
    verbose : bool, optional
        If True, print Slater integrals and Hilbert-space dimensions.
        Default is False.
    sparse_U : bool, optional
        If False (default), return dense rank-4 Coulomb tensors. If True,
        embed the impurity/core interactions directly as CSR matrices,
        avoiding dense rank-4 tensors over the full impurity and bath space.
        The flattened convention is ``row = lorb * norb + korb`` and
        ``col = jorb * norb + iorb``, where ``norb`` is ``ntot_v`` for the
        initial state and ``ntot`` for the intermediate state.
    tol : float, optional
        Keep Coulomb entries with absolute value greater than this threshold
        when ``sparse_U=True``. Unused for dense tensors. Default is ``1e-10``.

    Returns
    -------
    emat_i : numpy.ndarray, shape (ntot, ntot)
        Initial-state one-body matrix over impurity, bath, and core orbitals.
    umat_i : numpy.ndarray or scipy.sparse.csr_matrix
        Initial-state Coulomb interaction on impurity and core orbitals.
        Dense shape is ``(ntot,) * 4``; sparse shape is
        ``(ntot**2, ntot**2)``.
    basis_i : FockBasisSpec
        Initial-state basis specification with ``v_noccu`` electrons in
        the impurity and bath orbitals and ``c_norb`` core electrons.
    emat_n : numpy.ndarray, shape (ntot, ntot)
        Intermediate-state one-body matrix over impurity, bath, and core.
    umat_n : numpy.ndarray or scipy.sparse.csr_matrix
        Intermediate-state Coulomb interaction on impurity and core orbitals.
        Dense shape is ``(ntot,) * 4``; sparse shape is
        ``(ntot**2, ntot**2)``.
    basis_n : FockBasisSpec
        Intermediate-state basis specification with ``v_noccu + 1`` impurity
        and bath electrons and ``c_norb - 1`` core electrons.
    trans_mat : numpy.ndarray, shape (npol, ntot, ntot)
        Core-to-impurity absorption matrices, with nonzero entries only in
        the impurity-row/core-column block. ``npol`` is 3 for dipole and 5
        for quadrupole transitions; components are in the global frame.
    shift : float
        Filled core shell-level energy, ``c_norb * c_level``, excluding
        Coulomb contributions. Pass the first seven outputs to
        :func:`~edrixs.solvers.get_ops` and ``shift`` to dense, SciPy, or PETSc
        :func:`~edrixs.solvers.ed`. Integrals remain unshifted.
    """
    if verbose:
        print("edrixs >>> Setting up SIAM problem ...")

    v_name_options = ['s', 'p', 't2g', 'd', 'f']
    c_name_options = [
        's', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52',
        'f', 'f52', 'f72'
    ]

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

    if verbose:
        _print_slater_summary(slater_name, slater_i, slater_n)

    umat_tmp_i = get_umat_slater(v_name + c_name, *slater_i)
    umat_tmp_n = get_umat_slater(v_name + c_name, *slater_n)

    if sparse_U:
        umat_i = _embed_impurity_core_umat_sparse(
            umat_tmp_i, v_norb, c_norb, ntot_v, tol=tol
        )
        umat_n = _embed_impurity_core_umat_sparse(
            umat_tmp_n, v_norb, c_norb, ntot_v, tol=tol
        )
    else:
        umat_i = _embed_impurity_core_umat(
            umat_tmp_i, v_norb, c_norb, ntot_v
        )
        umat_n = _embed_impurity_core_umat(
            umat_tmp_n, v_norb, c_norb, ntot_v
        )

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)

    if siam_type == 1:
        if hopping is not None:
            emat_i[0:ntot_v, 0:ntot_v] += np.asarray(hopping)

        if hopping_n is not None:
            emat_n[0:ntot_v, 0:ntot_v] += np.asarray(hopping_n)
        elif hopping is not None:
            emat_n[0:ntot_v, 0:ntot_v] += np.asarray(hopping)

    elif siam_type == 0:
        if imp_mat is not None:
            emat_i[0:v_norb, 0:v_norb] += np.asarray(imp_mat)

        if imp_mat_n is not None:
            emat_n[0:v_norb, 0:v_norb] += np.asarray(imp_mat_n)
        elif imp_mat is not None:
            emat_n[0:v_norb, 0:v_norb] += np.asarray(imp_mat)

        if bath_level is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx = (i + 1) * v_norb + j
                    emat_i[idx, idx] += bath_level[i, j]

        if bath_level_n is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx = (i + 1) * v_norb + j
                    emat_n[idx, idx] += bath_level_n[i, j]
        elif bath_level is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx = (i + 1) * v_norb + j
                    emat_n[idx, idx] += bath_level[i, j]

        if hyb is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx1 = j
                    idx2 = (i + 1) * v_norb + j
                    emat_i[idx1, idx2] += hyb[i, j]
                    emat_i[idx2, idx1] += np.conj(hyb[i, j])

        if hyb_n is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx1 = j
                    idx2 = (i + 1) * v_norb + j
                    emat_n[idx1, idx2] += hyb_n[i, j]
                    emat_n[idx2, idx1] += np.conj(hyb_n[i, j])
        elif hyb is not None:
            for i in range(nbath):
                for j in range(v_norb):
                    idx1 = j
                    idx2 = (i + 1) * v_norb + j
                    emat_n[idx1, idx2] += hyb[i, j]
                    emat_n[idx2, idx1] += np.conj(hyb[i, j])

    else:
        raise Exception("Unknown siam_type: ", siam_type)

    if c_name in ['p', 'd', 'f']:
        emat_i[ntot_v:ntot, ntot_v:ntot] += atom_hsoc(c_name, c_soc)
        emat_n[ntot_v:ntot, ntot_v:ntot] += atom_hsoc(c_name, c_soc)

    # Static core-hole potential on the impurity in the intermediate state.
    emat_n[0:v_norb, 0:v_norb] -= np.eye(v_norb) * static_core_pot

    if trans_c2n is None:
        trans_c2n = np.eye(v_norb, dtype=complex)
    else:
        trans_c2n = np.asarray(trans_c2n)

    tmat = np.eye(ntot, dtype=complex)
    for i in range(nbath + 1):
        off = i * v_norb
        tmat[off:off + v_norb, off:off + v_norb] = np.conj(
            np.transpose(trans_c2n)
        )

    emat_i[:, :] = cb_op(emat_i, tmat)
    emat_n[:, :] = cb_op(emat_n, tmat)

    if ext_B is not None:
        zeeman = _valence_zeeman_matrix(v_name, v_orbl, ext_B, on_which)
        emat_i[0:v_norb, 0:v_norb] += zeeman
        emat_n[0:v_norb, 0:v_norb] += zeeman

    # Core levels in both sectors.
    emat_i[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level
    emat_n[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level

    basis_i = FockBasisSpec.from_args(ntot_v, v_noccu, c_norb, c_norb)
    basis_n = FockBasisSpec.from_args(ntot_v, v_noccu + 1, c_norb, c_norb - 1)

    if verbose:
        print("edrixs >>> Dimension of the initial Hamiltonian: ", len(basis_i))
        print("edrixs >>> Dimension of the intermediate Hamiltonian: ", len(basis_n))

    tmp_g = _rotated_transition_blocks(v_name + c_name, loc_axis)
    trans_mat = np.zeros((tmp_g.shape[0], ntot, ntot), dtype=complex)
    trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

    if verbose:
        print("edrixs >>> SIAM setup Done !")

    shift = float(c_norb * c_level)
    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, shift


def model_siam_2d1p(
    *, slater, impurity_levels, bath_levels, hyb, Delta, nd,
    v_soc=None, c_soc=0, om_shift=0, ext_B=None, on_which='spin',
    loc_axis=None, verbose=False, sparse_U=False, tol=1E-10
):
    """
    Set up the commonly used d-impurity, one-bath, p-core charge-transfer
    model. This is a convenience wrapper around :func:`model_siam`.
    All energies, Slater integrals, hoppings, SOC constants, and magnetic-field
    components are in eV. The bath has ten spin-orbitals and no Coulomb interaction.

    Parameters
    ----------
    slater : tuple of array_like
        Initial ``[F0_dd, F2_dd, F4_dd]`` and intermediate
        ``[F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp]``.
        Each state's average Coulomb interactions are recovered by subtracting
        the corresponding ``get_F0`` contribution from its F0 integrals.
    impurity_levels, bath_levels : array_like, shape (5,)
        Real orbital offsets in the order ``(dz2, dzx, dzy, dx2-y2, dxy)``.
        The mean of each five-orbital array is subtracted before adding the
        charge-transfer-derived shell center. Each level applies to both spins;
        the input arrays are not modified.
    hyb : array_like, shape (5,)
        Orbital-diagonal impurity-to-bath hoppings in the same order. Complex
        values are allowed. Both spins and both states use the same hoppings.
    Delta : float
        Charge-transfer energy. This is defined in the atomic limit
        (without crystal field, hybridization, or multipolar Coulomb
        interactions) as the energy cost of transferring
        one electron from the bath to the impurity.
        Used for both states in :func:`~edrixs.utils.CT_imp_bath` and
        :func:`~edrixs.utils.CT_imp_bath_core_hole`.
    nd : int
        Nominal impurity occupancy, from 0 through 9. The nominal bath is full;
        the initial total valence occupancy is ``nd + 10``.
    v_soc : pair of float, optional
        Initial and intermediate impurity SOC constants. None means zero.
    c_soc : float, optional
        Core SOC constant. Default zero.
    om_shift : float, optional
        Overall XAS energy alignment. The core energy passed to model_siam is
        ``-om_shift - 5*E_p``, where E_p comes from the core-hole CT helper.
        Increasing om_shift raises transition energies by the same amount.
        Default zero; the filled-core reference convention is preserved.
    ext_B, on_which, loc_axis, verbose, sparse_U, tol
        Passed through to :func:`model_siam` with the same meaning and defaults.
        The magnetic field acts on the impurity only.

    Returns
    -------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, shift
        The eight model_siam outputs. Pass the first seven to
        :func:`~edrixs.solvers.get_ops` and the final ``shift`` to dense, SciPy, or
        PETSc :func:`~edrixs.solvers.ed`. Here ``shift = 6 * (-om_shift - 5*E_p)``;
        it excludes Coulomb contributions and does not alter the integrals. One-body
        matrices both have shape (26, 26) and use complex spherical
        harmonics. Basis metadata fixes initial valence occupancy to nd+10,
        intermediate valence occupancy to nd+11, and core occupancies to six
        and five, respectively.
        Coulomb tensors are dense unless sparse_U=True. No diagonalization
        or spectrum calculation is performed.
    """
    slater_i, slater_n = slater
    if v_soc is None:
        v_soc = (0, 0)

    # Center the orbital splittings before applying the CT-derived energies.
    impurity_levels = np.asarray(impurity_levels) - np.mean(impurity_levels)
    bath_levels = np.asarray(bath_levels) - np.mean(bath_levels)
    U_dd_i = slater_i[0] - get_F0('d', *slater_i[1:3])
    U_dd_n = slater_n[0] - get_F0('d', *slater_n[1:3])
    U_dp = slater_n[3] - get_F0('dp', *slater_n[5:7])
    E_d, E_L = CT_imp_bath(U_dd_i, Delta, nd)
    E_dc, E_Lc, E_p = CT_imp_bath_core_hole(U_dd_n, U_dp, Delta, nd)

    trans_c2n = tmat_c2r('d', True)
    imp_mat = (np.diag(np.repeat(impurity_levels, 2) + E_d)
               + cb_op(atom_hsoc('d', v_soc[0]), trans_c2n))
    imp_mat_n = (np.diag(np.repeat(impurity_levels, 2) + E_dc)
                 + cb_op(atom_hsoc('d', v_soc[1]), trans_c2n))
    bath_offsets = np.repeat(bath_levels, 2)[None, :]
    return model_siam(
        ('d', 'p'), 1, siam_type=0, v_noccu=nd + 10,
        c_level=-om_shift - 5 * E_p, c_soc=c_soc, trans_c2n=trans_c2n,
        imp_mat=imp_mat, imp_mat_n=imp_mat_n,
        bath_level=bath_offsets + E_L, bath_level_n=bath_offsets + E_Lc,
        hyb=np.repeat(hyb, 2)[None, :], slater=(slater_i, slater_n),
        ext_B=ext_B, on_which=on_which, loc_axis=loc_axis,
        verbose=verbose, sparse_U=sparse_U, tol=tol,
    )
