"""Backend-independent physical model construction for EDRIXS.

The model functions in this module describe orbital-space Hamiltonians,
interactions, Fock bases, and photon-transition matrices without constructing
backend-owned many-body operators.
"""

from __future__ import annotations

import numpy as np

from .angular_momentum import get_lx, get_ly, get_lz, get_sx, get_sy, get_sz
from .angular_momentum import get_wigner_dmat, rmat_to_euler
from .basis_transform import cb_op, tmat_r2c
from .coulomb_utensor import get_umat_slater, get_umat_slater_3shells
from .fock_basis import get_fock_basis_int
from .iostream import write_emat, write_umat
from .photon_transition import get_trans_oper
from .soc import atom_hsoc
from .utils import info_atomic_shell, slater_integrals_name
from ._solvers_helpers import (
    _embed_impurity_core_umat,
    _embed_impurity_core_umat_sparse,
    _rotated_transition_blocks,
    _umat_dense_to_sparse,
    _valence_zeeman_matrix,
)

__all__ = ['model_1v1c', 'model_2v1c', 'model_siam']


def model_1v1c(shell_name, *, shell_level=None, v_soc=None, c_soc=0,
               v_noccu=1, slater=None, ext_B=None, on_which='spin',
               v_cfmat=None, v_othermat=None, loc_axis=None, verbose=0,
               sparse_U=False, tol=1E-10):
    """
    Set up orbital-space data and Fock bases for a 1v1c problem.

    This routine defines the physical one-valence-shell/one-core-shell problem
    independently of the numerical backend. It constructs the one-body orbital
    matrices, Coulomb tensors, Fock bases, and orbital-space transition
    matrices, but it does not build many-body Hamiltonians and does not
    diagonalize anything.

    Parameters
    ----------
    shell_name : tuple of str
        Names of the valence and core shells.

    shell_level, v_soc, c_soc, v_noccu, slater, ext_B, on_which,
    v_cfmat, v_othermat, loc_axis, verbose
        Same physical model parameters as ed_1v1c_py.

    sparse_U : bool, optional
        If False, return dense rank-4 Coulomb tensors. If True, return each
        Coulomb tensor as a sparse matrix with shape (ntot*ntot, ntot*ntot),
        using the flattened convention
        row = lorb * ntot + korb and col = jorb * ntot + iorb.

    tol : float, optional
        Threshold used when converting dense Coulomb tensors to sparse format.

    Returns
    -------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat
        Backend-independent problem definition. trans_mat has shape
        (npol, ntot, ntot).
    """
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

    print()
    print("    Summary of Slater integrals:")
    print("    ------------------------------")
    print("    Terms,   Initial Hamiltonian,  Intermediate Hamiltonian")
    for i in range(nslat):
        print("    ", slater_name[i], ":  {:20.10f}{:20.10f}".format(
            slater_i[i], slater_n[i]
        ))
    print()

    case = v_name + c_name
    umat_i = get_umat_slater(case, *slater_i)
    umat_n = get_umat_slater(case, *slater_n)

    if verbose > 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

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
        emat_i[v_norb:ntot, v_norb:ntot] += np.eye(c_norb) * shell_level[1]

        emat_n[0:v_norb, 0:v_norb] += np.eye(v_norb) * shell_level[0]
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

    if verbose > 0:
        write_emat(emat_i, 'hopping_i.in')
        write_emat(emat_n, 'hopping_n.in')

    # Fock bases.
    basis_i = get_fock_basis_int(v_norb, v_noccu, c_norb, c_norb)
    basis_n = get_fock_basis_int(v_norb, v_noccu + 1, c_norb, c_norb - 1)

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

    print("edrixs >>> 1v1c setup Done !")

    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat


def model_2v1c(
    shell_name, *, shell_level=None,
    v1_soc=None, v2_soc=None, c_soc=0, v_tot_noccu=1, slater=None,
    v1_ext_B=None, v2_ext_B=None, v1_on_which='spin',
    v2_on_which='spin', v1_cfmat=None, v2_cfmat=None,
    v1_othermat=None, v2_othermat=None, hopping_v1v2=None,
    trans_to_which=1, loc_axis=None, verbose=0, sparse_U=False, tol=1E-10
):
    """
    Set up orbital-space data and Fock bases for a 2v1c problem.

    This is the backend-neutral setup analogue of the 2-valence-shell,
    1-core-shell Fortran ED/RIXS input construction.  It does not build
    many-body Hamiltonians and does not diagonalize anything. If sparse_U=True,
    the returned Coulomb tensors are sparse flattened matrices rather than
    dense rank-4 arrays.

    Returns
    -------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat
        These can be passed directly to get_ops(..., backend='scipy') or
        get_ops(..., backend='dense').
    """
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

    print()
    print("    Summary of Slater integrals:")
    print("    ------------------------------")
    print("    Terms,   Initial Hamiltonian,  Intermediate Hamiltonian")
    for i in range(nslat):
        print("    ", slater_name[i], ":  {:20.10f}{:20.10f}".format(
            slater_i[i], slater_n[i]
        ))
    print()

    umat_i = get_umat_slater_3shells(
        (v1_name, v2_name, c_name), *slater_i
    )
    umat_n = get_umat_slater_3shells(
        (v1_name, v2_name, c_name), *slater_n
    )

    if verbose > 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

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

    # Shell levels.  Since the setup/get_ops route keeps the core orbitals in the
    # Fock basis, the filled core contribution is represented explicitly.
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

    if verbose > 0:
        write_emat(emat_i, 'hopping_i.in')
        write_emat(emat_n, 'hopping_n.in')

    basis_i = get_fock_basis_int(
        v1v2_norb, v_tot_noccu, c_norb, c_norb
    )
    basis_n = get_fock_basis_int(
        v1v2_norb, v_tot_noccu + 1, c_norb, c_norb - 1
    )

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

    print("edrixs >>> 2v1c setup Done !")

    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat


def model_siam(
    shell_name, nbath, *, siam_type=0, v_noccu=1, static_core_pot=0,
    c_level=0, c_soc=0, trans_c2n=None, imp_mat=None, imp_mat_n=None,
    bath_level=None, bath_level_n=None, hyb=None, hyb_n=None,
    hopping=None, hopping_n=None, slater=None, ext_B=None,
    on_which='spin', loc_axis=None, verbose=0, sparse_U=False, tol=1E-10
):
    """
    Set up orbital-space data and Fock bases for a SIAM problem.

    This is the backend-neutral setup analogue of ed_siam_fort. It does not
    search over occupancies, does not build many-body Hamiltonians, and does not
    diagonalize anything. The occupancy is the supplied v_noccu. If
    sparse_U=True, the impurity+core Coulomb tensor is embedded directly into
    the full SIAM orbital space as a sparse flattened matrix.

    Returns
    -------
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat
        These can be passed directly to get_ops(..., backend='scipy') or
        get_ops(..., backend='dense').
    """
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

    print()
    print("    Summary of Slater integrals:")
    print("    ------------------------------")
    print("    Terms,   Initial Hamiltonian,  Intermediate Hamiltonian")
    for i in range(nslat):
        print("    ", slater_name[i], ":  {:20.10f}{:20.10f}".format(
            slater_i[i], slater_n[i]
        ))
    print()

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

    if verbose > 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

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

    # Since the setup/get_ops route keeps the core orbitals in the Fock basis, the
    # filled-core contribution is represented explicitly.
    emat_i[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level
    emat_n[ntot_v:ntot, ntot_v:ntot] += np.eye(c_norb) * c_level

    if verbose > 0:
        write_emat(emat_i, 'hopping_i.in')
        write_emat(emat_n, 'hopping_n.in')

    basis_i = get_fock_basis_int(ntot_v, v_noccu, c_norb, c_norb)
    basis_n = get_fock_basis_int(ntot_v, v_noccu + 1, c_norb, c_norb - 1)

    print("edrixs >>> Dimension of the initial Hamiltonian: ", len(basis_i))
    print("edrixs >>> Dimension of the intermediate Hamiltonian: ", len(basis_n))

    tmp_g = _rotated_transition_blocks(v_name + c_name, loc_axis)
    trans_mat = np.zeros((tmp_g.shape[0], ntot, ntot), dtype=complex)
    trans_mat[:, 0:v_norb, ntot_v:ntot] = tmp_g

    print("edrixs >>> SIAM setup Done !")

    return emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat
