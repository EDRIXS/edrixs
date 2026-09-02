"""Private implementation helpers for :mod:`edrixs.solvers`."""

import numpy as np
import scipy.sparse as sp

from .iostream import write_emat, write_umat, write_config, read_poles_from_file
from .angular_momentum import (
    get_sx, get_sy, get_sz, get_lx, get_ly, get_lz, rmat_to_euler, get_wigner_dmat
)
from .photon_transition import (
    get_trans_oper, quadrupole_polvec, dipole_polvec_xas, dipole_polvec_rixs, unit_wavevector
)
from .coulomb_utensor import get_umat_slater, get_umat_slater_3shells
from .fock_basis import write_fock_dec_by_N
from .basis_transform import tmat_r2c
from .utils import info_atomic_shell, slater_integrals_name
from .plot_spectrum import get_spectra_from_poles, merge_pole_dicts
from .soc import atom_hsoc


def _infer_backend(*operators):
    """Infer which backend owns every supplied operator."""
    from .petsc_backend import petsc_backend
    from .scipy_backend import scipy_backend
    from .fortran_backend import fortran_backend

    owned_by_scipy = bool(operators) and all(
        scipy_backend.owns_operator_scipy(operator) for operator in operators
    )
    owned_by_petsc = bool(operators) and all(
        petsc_backend.owns_operator_petsc(operator) for operator in operators
    )
    owned_by_fortran = bool(operators) and all(
        fortran_backend.owns_operator_fortran(operator) for operator in operators
    )

    match owned_by_scipy, owned_by_petsc, owned_by_fortran:
        case True, False, False:
            return 'scipy'
        case False, True, False:
            return 'petsc'
        case False, False, True:
            return 'fortran'
        case False, False, False:
            raise TypeError(
                "Could not infer a backend from the supplied operators; "
                "pass backend='scipy', backend='petsc', or backend='fortran' explicitly"
            )
        case _:
            raise TypeError(
                "Backend inference is ambiguous; pass backend explicitly"
            )


def _rotated_transition_blocks(case, loc_axis=None):
    """
    Build transition blocks in global coordinates for a shell-transition case.
    """
    if loc_axis is None:
        loc_axis = np.eye(3)
    else:
        loc_axis = np.asarray(loc_axis)

    tmp = get_trans_oper(case)
    npol, n, m = tmp.shape
    tmp_g = np.zeros((npol, n, m), dtype=complex)

    if npol == 3:
        for i in range(3):
            for j in range(3):
                tmp_g[i] += loc_axis[i, j] * tmp[j]
    elif npol == 5:
        alpha, beta, gamma = rmat_to_euler(loc_axis)
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

    return tmp_g


def _valence_zeeman_matrix(shell_name, orbl, ext_B, on_which):
    """
    Build a Zeeman matrix for a valence shell.
    """
    if shell_name == 't2g':
        lx, ly, lz = get_lx(1, True), get_ly(1, True), get_lz(1, True)
        sx, sy, sz = get_sx(1), get_sy(1), get_sz(1)
        lx, ly, lz = -lx, -ly, -lz
    else:
        lx, ly, lz = get_lx(orbl, True), get_ly(orbl, True), get_lz(orbl, True)
        sx, sy, sz = get_sx(orbl), get_sy(orbl), get_sz(orbl)

    if on_which.strip() == 'spin':
        return ext_B[0] * (2 * sx) + ext_B[1] * (2 * sy) + ext_B[2] * (2 * sz)

    if on_which.strip() == 'orbital':
        return ext_B[0] * lx + ext_B[1] * ly + ext_B[2] * lz

    if on_which.strip() == 'both':
        return (
            ext_B[0] * (lx + 2 * sx)
            + ext_B[1] * (ly + 2 * sy)
            + ext_B[2] * (lz + 2 * sz)
        )

    raise Exception("Unknown value of on_which", on_which)


def _umat_dense_to_sparse(umat, tol=1E-10):
    """
    Convert a dense rank-4 Coulomb tensor to a sparse flattened matrix.

    The flattening convention is

        row = lorb * norbs + korb
        col = jorb * norbs + iorb

    corresponding to tensor entry umat[lorb, korb, jorb, iorb].
    """
    umat = np.asarray(umat)

    if umat.ndim != 4:
        raise ValueError("dense umat must be a rank-4 tensor")

    norbs = umat.shape[0]
    if umat.shape != (norbs, norbs, norbs, norbs):
        raise ValueError("dense umat must have shape (n, n, n, n)")

    lorb, korb, jorb, iorb = np.nonzero(np.abs(umat) > tol)

    rows = lorb * norbs + korb
    cols = jorb * norbs + iorb
    data = umat[lorb, korb, jorb, iorb]

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(norbs * norbs, norbs * norbs),
        dtype=np.complex128
    ).tocsr()


def _embed_impurity_core_umat_sparse(umat_tmp, v_norb, c_norb, ntot_v,
                                     tol=1E-10):
    """
    Embed an impurity+core Coulomb tensor into full SIAM space sparsely.

    The input umat_tmp is the dense tensor for the compact impurity+core
    orbital space of size v_norb + c_norb. The returned object is a sparse
    flattened matrix for the full SIAM space of size ntot_v + c_norb.

    No dense full-space rank-4 tensor is allocated.
    """
    umat_tmp = np.asarray(umat_tmp)

    nsmall = v_norb + c_norb
    if umat_tmp.shape != (nsmall, nsmall, nsmall, nsmall):
        raise ValueError(
            "umat_tmp has shape {}, expected {}".format(
                umat_tmp.shape, (nsmall, nsmall, nsmall, nsmall)
            )
        )

    ntot = ntot_v + c_norb
    index_map = np.array(
        list(range(v_norb)) + [ntot_v + i for i in range(c_norb)],
        dtype=int
    )

    a, b, c, d = np.nonzero(np.abs(umat_tmp) > tol)

    lorb = index_map[a]
    korb = index_map[b]
    jorb = index_map[c]
    iorb = index_map[d]

    rows = lorb * ntot + korb
    cols = jorb * ntot + iorb
    data = umat_tmp[a, b, c, d]

    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(ntot * ntot, ntot * ntot),
        dtype=np.complex128
    ).tocsr()


def _embed_impurity_core_umat(umat_tmp, v_norb, c_norb, ntot_v):
    """
    Embed an impurity+core Coulomb tensor into the full SIAM orbital space.

    The bath orbitals do not carry local Coulomb interactions in this SIAM
    representation, so only impurity and core indices are populated.
    """
    ntot = ntot_v + c_norb
    umat = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)

    idx = list(range(v_norb)) + [ntot_v + i for i in range(c_norb)]

    for i in range(v_norb + c_norb):
        for j in range(v_norb + c_norb):
            for k in range(v_norb + c_norb):
                for m in range(v_norb + c_norb):
                    umat[idx[i], idx[j], idx[k], idx[m]] = umat_tmp[i, j, k, m]

    return umat


def _expand_broadening(gamma, n, name):
    """Return gamma as a length-n float array."""
    out = np.empty(n, dtype=float)
    if np.isscalar(gamma):
        out[:] = gamma
    else:
        gamma = np.asarray(gamma, dtype=float)
        if gamma.shape != (n,):
            raise ValueError("{} must be scalar or have shape ({},)".format(name, n))
        out[:] = gamma
    return out


def _ed_1or2_valence_1core(
        comm, shell_name, *, shell_level=None,
        v1_soc=None, v2_soc=None, c_soc=0, v_tot_noccu=1, slater=None,
        v1_ext_B=None, v2_ext_B=None, v1_on_which='spin', v2_on_which='spin',
        v1_cfmat=None, v2_cfmat=None, v1_othermat=None, v2_othermat=None,
        hopping_v1v2=None, do_ed=True, ed_solver=2, neval=1, nvector=1, ncv=3,
        idump=False, maxiter=500, eigval_tol=1e-8, min_ndim=1000
        ):
    from .fedrixs import ed_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()
    if rank == 0:
        print("edrixs >>> Running ED ...", flush=True)
    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()
    info_shell = info_atomic_shell()

    # Quantum numbers of angular momentum
    v1_orbl = info_shell[v1_name][0]
    if v2_name != 'empty':
        v2_orbl = info_shell[v2_name][0]
    else:
        v2_orbl = -1

    # number of orbitals with spin
    v1_norb = info_shell[v1_name][1]
    if v2_name != 'empty':
        v2_norb = info_shell[v2_name][1]
    else:
        v2_norb = 0
    c_norb = info_shell[c_name][1]

    # total number of orbitals
    ntot = v1_norb + v2_norb + c_norb
    v1v2_norb = v1_norb + v2_norb

    # Coulomb interaction
    if v2_name == 'empty':
        slater_name = slater_integrals_name((v1_name, c_name), ('v', 'c'))
    else:
        slater_name = slater_integrals_name((v1_name, v2_name, c_name), ('v1', 'v2', 'c1'))
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

    if v2_name == 'empty':
        umat_i = get_umat_slater(v1_name + c_name, *slater_i)
        umat_n = get_umat_slater(v1_name + c_name, *slater_n)
    else:
        umat_i = get_umat_slater_3shells((v1_name, v2_name, c_name), *slater_i)
        umat_n = get_umat_slater_3shells((v1_name, v2_name, c_name), *slater_n)

    if rank == 0:
        write_umat(umat_i, 'coulomb_i.in')
        write_umat(umat_n, 'coulomb_n.in')

    emat_i = np.zeros((ntot, ntot), dtype=complex)
    emat_n = np.zeros((ntot, ntot), dtype=complex)
    # SOC
    if v1_soc is not None and v1_name in ['p', 'd', 't2g', 'f']:
        emat_i[0:v1_norb, 0:v1_norb] += atom_hsoc(v1_name, v1_soc[0])
        emat_n[0:v1_norb, 0:v1_norb] += atom_hsoc(v1_name, v1_soc[1])

    if v2_soc is not None and v2_name in ['p', 'd', 't2g', 'f']:
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += atom_hsoc(v2_name, v2_soc[0])
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += atom_hsoc(v2_name, v2_soc[1])

    if c_name in ['p', 'd', 'f']:
        emat_n[v1v2_norb:ntot, v1v2_norb:ntot] += atom_hsoc(c_name, c_soc)

    # crystal field
    if v1_cfmat is not None:
        emat_i[0:v1_norb, 0:v1_norb] += np.array(v1_cfmat)
        emat_n[0:v1_norb, 0:v1_norb] += np.array(v1_cfmat)

    if v2_cfmat is not None and v2_name != 'empty':
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.array(v2_cfmat)
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.array(v2_cfmat)

    # other mat
    if v1_othermat is not None:
        emat_i[0:v1_norb, 0:v1_norb] += np.array(v1_othermat)
        emat_n[0:v1_norb, 0:v1_norb] += np.array(v1_othermat)

    if v2_othermat is not None and v2_name != 'empty':
        emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.array(v2_othermat)
        emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.array(v2_othermat)

    # energy of shell
    if shell_level is not None:
        eval_shift = shell_level[2] * c_norb / v_tot_noccu
        emat_i[0:v1_norb, 0:v1_norb] += np.eye(v1_norb) * shell_level[0]
        emat_i[0:v1_norb, 0:v1_norb] += np.eye(v1_norb) * eval_shift
        emat_n[0:v1_norb, 0:v1_norb] += np.eye(v1_norb) * shell_level[0]
        emat_n[v1v2_norb:ntot, v1v2_norb:ntot] += np.eye(c_norb) * shell_level[2]
        if v2_name != 'empty':
            emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.eye(v2_norb) * shell_level[1]
            emat_i[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.eye(v2_norb) * eval_shift
            emat_n[v1_norb:v1v2_norb, v1_norb:v1v2_norb] += np.eye(v2_norb) * shell_level[1]

    # external magnetic field
    for name, l, ext_B, which, i1, i2 in [
        (v1_name, v1_orbl, v1_ext_B, v1_on_which, 0, v1_norb),
        (v2_name, v2_orbl, v2_ext_B, v2_on_which, v1_norb, v1v2_norb)
    ]:
        if name == 'empty':
            continue
        if name == 't2g':
            lx, ly, lz = get_lx(1, True), get_ly(1, True), get_lz(1, True)
            sx, sy, sz = get_sx(1), get_sy(1), get_sz(1)
            lx, ly, lz = -lx, -ly, -lz
        else:
            lx, ly, lz = get_lx(l, True), get_ly(l, True), get_lz(l, True)
            sx, sy, sz = get_sx(l), get_sy(l), get_sz(l)
        if ext_B is not None:
            if which.strip() == 'spin':
                zeeman = ext_B[0] * (2 * sx) + ext_B[1] * (2 * sy) + ext_B[2] * (2 * sz)
            elif which.strip() == 'orbital':
                zeeman = ext_B[0] * lx + ext_B[1] * ly + ext_B[2] * lz
            elif which.strip() == 'both':
                zeeman = (ext_B[0] * (lx + 2 * sx) +
                          ext_B[1] * (ly + 2 * sy) +
                          ext_B[2] * (lz + 2 * sz))
            else:
                raise Exception("Unknown value of zeeman_on_which", which)
            emat_i[i1:i2, i1:i2] += zeeman
            emat_n[i1:i2, i1:i2] += zeeman

    # hopping between the two valence shells
    if hopping_v1v2 is not None and v2_name != 'empty':
        emat_i[0:v1_norb, v1_norb:v1v2_norb] += np.array(hopping_v1v2)
        emat_i[v1_norb:v1v2_norb, 0:v1_norb] += np.conj(np.transpose(hopping_v1v2))
        emat_n[0:v1_norb, v1_norb:v1v2_norb] += np.array(hopping_v1v2)
        emat_n[v1_norb:v1v2_norb, 0:v1_norb] += np.conj(np.transpose(hopping_v1v2))

    if rank == 0:
        write_emat(emat_i, 'hopping_i.in')
        write_emat(emat_n, 'hopping_n.in')
        write_config(
            './', ed_solver, v1v2_norb, c_norb, neval, nvector, ncv, idump,
            maxiter=maxiter, min_ndim=min_ndim, eigval_tol=eigval_tol
        )
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu, "fock_i.in")

    if do_ed:
        # now, call ed solver
        comm.Barrier()
        ed_fsolver(fcomm, rank, size)
        comm.Barrier()

        # read eigvals.dat and denmat.dat
        data = np.loadtxt('eigvals.dat', ndmin=2)
        eval_i = np.zeros(neval, dtype=float)
        eval_i[0:neval] = data[0:neval, 1]
        data = np.loadtxt('denmat.dat', ndmin=2)
        tmp = (nvector, v1v2_norb, v1v2_norb)
        denmat = data[:, 3].reshape(tmp) + 1j * data[:, 4].reshape(tmp)

        return eval_i, denmat
    else:
        return None, None


def _xas_1or2_valence_1core(
        comm, shell_name, ominc, *, gamma_c=0.1,
        v_tot_noccu=1, trans_to_which=1, thin=1.0, phi=0,
        pol_type=None, num_gs=1, nkryl=200, temperature=1.0,
        loc_axis=None, scatter_axis=None
        ):
    from .fedrixs import xas_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()

    info_shell = info_atomic_shell()
    v1_norb = info_shell[v1_name][1]
    if v2_name != 'empty':
        v2_norb = info_shell[v2_name][1]
    else:
        v2_norb = 0

    c_norb = info_shell[c_name][1]
    ntot = v1_norb + v2_norb + c_norb
    v1v2_norb = v1_norb + v2_norb
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
        write_config(num_val_orbs=v1v2_norb, num_core_orbs=c_norb,
                     num_gs=num_gs, nkryl=nkryl)
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu, "fock_i.in")
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu + 1, "fock_n.in")

    # Build transition operators in local-xyz axis
    if trans_to_which == 1:
        case = v1_name + c_name
    elif trans_to_which == 2 and v2_name != 'empty':
        case = v2_name + c_name
    else:
        raise Exception('Unkonwn trans_to_which: ', trans_to_which)
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
    if trans_to_which == 1:
        trans_mat[:, 0:v1_norb, v1v2_norb:ntot] = tmp_g
    else:
        trans_mat[:, v1_norb:v1v2_norb, v1v2_norb:ntot] = tmp_g

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


def _rixs_1or2_valence_1core(
        comm, shell_name, ominc, eloss, *, gamma_c=0.1, gamma_f=0.1,
        v_tot_noccu=1, trans_to_which=1, thin=1.0, thout=1.0, phi=0,
        pol_type=None, num_gs=1, nkryl=200, linsys_max=500, linsys_tol=1e-8,
        temperature=1.0, loc_axis=None, scatter_axis=None
        ):
    from .fedrixs import rixs_fsolver

    rank = comm.Get_rank()
    size = comm.Get_size()
    fcomm = comm.py2f()

    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    c_name = shell_name[2].strip()

    info_shell = info_atomic_shell()

    v1_norb = info_shell[v1_name][1]
    if v2_name != 'empty':
        v2_norb = info_shell[v2_name][1]
    else:
        v2_norb = 0
    c_norb = info_shell[c_name][1]
    ntot = v1_norb + v2_norb + c_norb
    v1v2_norb = v1_norb + v2_norb
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
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu, "fock_i.in")
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu + 1, "fock_n.in")
        write_fock_dec_by_N(v1v2_norb, v_tot_noccu, "fock_f.in")

        # Build transition operators in local-xyz axis
        if trans_to_which == 1:
            case = v1_name + c_name
        elif trans_to_which == 2:
            case = v2_name + c_name
        else:
            raise Exception('Unkonwn trans_to_which: ', trans_to_which)
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
        if trans_to_which == 1:
            trans_mat[:, 0:v1_norb, v1v2_norb:ntot] = tmp_g
        else:
            trans_mat[:, v1_norb:v1v2_norb, v1v2_norb:ntot] = tmp_g

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
                num_val_orbs=v1v2_norb, num_core_orbs=c_norb,
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
