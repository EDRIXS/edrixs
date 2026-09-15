import numpy as np
import pandas as pd
import edrixs
from itertools import permutations
import h5py

from math import comb


# Conversion to Sunny

A_RYABOV = {(1, 0): 1, (1, 1): 1,
            (2, 0): 3, (2, 1): 1, (2, 2): 1,
            (3, 0): 5, (3, 1): 5, (3, 2): 1, (3, 3): 1,
            (4, 0): 35, (4, 1): 7, (4, 2): 7, (4, 3): 1, (4, 4): 1,
            (5, 0): 63, (5, 1): 21, (5, 2): 3, (5, 3): 9, (5, 4): 1, (5, 5): 1,
            (6, 0): 231, (6, 1): 33, (6, 2): 33, (6, 3): 11, (6, 4): 11, (6, 5): 1, (6, 6): 1}


def stevens_to_sunny(k, q):
    """
    Parameters
    ----------
    k : rank
    q : projection

    Returns
    -------
    c_kq(k,q)/A_ryabov(k,q)

    """
    q = abs(q)
    return np.sqrt((2 - (q == 0)) * comb(2 * k, k - q)) / A_RYABOV[(k, q)]

# Spin, quadrupoles, octupoles, and point groups


def spin_ops(S):
    """
    Parameters
    ----------
    S : spin value

    Returns
    -------
    Sx, Sy, Sz : spin matrices for spin-S

    Sp : Splus (Raising operator)
    Sm : Sminus (Lowering operator)

    """
    d = int(round(2 * S + 1))
    m = np.arange(S, -S - 1, -1)
    Sz = np.diag(m).astype(complex)
    Sp = np.zeros((d, d), complex)
    for i in range(d - 1):
        Sp[i, i + 1] = np.sqrt(S * (S + 1) - m[i + 1] * (m[i + 1] + 1))
    Sm = Sp.conj().T
    return 0.5 * (Sp + Sm), (Sp - Sm) / (2j), Sz, Sp, Sm          # Sx, Sy, Sz, S+, S-


def dipole_irrep_map(group='spherical'):
    """
    Parameters
    ----------
    group

    description
    -----------
        spherical, Oh                         :  t1g = (x, y, z)                          all three degenerate
        D6h, D4h, D3d, D2d, C4v, C3v, D4, C4h :
            a2g = (z,)   eg = (x, y)   uniaxial: parallel / in-plane
        D2h, D2                               :  b1g = (z,)   b2g = (y,)  b3g = (x,)      non-degenerate
        C1                                    :  x   = (x,)   y   = (y,)  z   = (z,)

    Returns
    -------
    dict = {irrep: (component indices into sigma1)}, with sigma1 ordered
    (x, y, z) as fundamental_spectra returns it.

    """
    g = (group or 'spherical').strip()
    if g in ('spherical', 'SO3', 'sph', 'Oh', 'O', 'Td'):
        return {'t1g': (0, 1, 2)}
    if g in ('D6h', 'D4h', 'D3d', 'D2d', 'C4v', 'C3v', 'D4', 'C4h'):
        return {'a2g': (2,), 'eg': (0, 1)}
    if g in ('D2h', 'D2'):
        return {'b1g': (2,), 'b2g': (1,), 'b3g': (0,)}
    if g in ('C1',):
        return {'x': (0,), 'y': (1,), 'z': (2,)}


def quadrupole_ops(S_low):
    """
    Parameters
    ----------
    S_low : spin operators

    Description
    -----------
    Q^{mu,nu} = (1/2){S_mu, S_nu} - (1/3) delta_{mu,nu} S^2

    Returns
    -------
    Q
    """
    Sv = np.asarray(S_low)
    S2 = sum(Sv[a] @ Sv[a] for a in range(3))
    d = Sv.shape[1]
    Q = np.zeros((3, 3, d, d), dtype=complex)
    for m in range(3):
        for n in range(3):
            Q[m, n] = 0.5 * (Sv[m] @ Sv[n] + Sv[n] @ Sv[m]) - (m == n) * S2 / 3.0
    return Q


def quadrupole_pol_basis(group='spherical'):
    """
    Parameters
    ----------
    group

    Description
    -----------
    The K=2 polarisation channel is 5-dimensional.  In a crystal it splits, and
    each piece carries its own fundamental spectrum:

        spherical : {K2}                                      one sigma^(2)
        Oh, O, Td        : eg (2) + t2g (3)                   two:  sigma^(2)_eg, sigma^(2)_t2g
        D4h, C4v, D4, C4h       : a1g + b1g + b2g + eg (2)    four
        D2h, D2 : a1g + a1g + b1g + b2g + b3g                 five

    Returns
    -------
    dict = {irrep: [E_1, E_2, ...]} with each E orthonormal in the Frobenius inner product.

    """
    s6, s2 = np.sqrt(6), np.sqrt(2)
    E_3z2 = np.diag([1., 1, -2]) / s6                      # 3z^2 - r^2
    E_x2y2 = np.diag([1., -1, 0]) / s2                     # x^2 - y^2
    off = {}
    for nm, (a, b) in {'xy': (0, 1), 'zx': (0, 2), 'yz': (1, 2)}.items():
        M = np.zeros((3, 3))
        M[a, b] = M[b, a] = 1
        off[nm] = M / s2

    g = (group or 'spherical').strip()
    if g in ('spherical', 'SO3', 'sph', 'C1'):
        return {'K2': [E_3z2, E_x2y2, off['xy'], off['zx'], off['yz']]}
    if g in ('Oh', 'O', 'Td'):
        return {'eg': [E_3z2, E_x2y2], 't2g': [off['xy'], off['zx'], off['yz']]}
    if g in ('D4h', 'C4v', 'D4', 'C4h'):
        return {'a1g': [E_3z2], 'b1g': [E_x2y2], 'b2g': [off['xy']], 'eg': [off['zx'], off['yz']]}
    if g in ('D2h', 'D2'):
        return {'a1g': [E_3z2], "a1g'": [E_x2y2], 'b1g': [off['xy']], 'b2g': [off['zx']], 'b3g': [off['yz']]}


def octupole_ops(S_low):
    """
    Parameters
    ----------
    S_low : spin operators

    Description
    -----------
    Rank-3 (octupole) Cartesian tensor operator, the rank-3 companion of quadrupole_ops():

        O^{mu,nu,rho} = sym(S_mu S_nu S_rho) - (1/5)[ d_{mu,nu} V_rho + d_{nu,rho} V_mu + d_{rho,mu} V_nu ],
        V_rho = sum_mu sym(S_mu S_mu S_rho),

    where sym(...) is the average over all 3! orderings.  The O-tensor is traceless in every pair of indices.

    Because O^{mu,nu,rho} is already fully symmetric and traceless, each
    tesseral operator is just the corresponding monomial combination of its
    Cartesian components -- no further projection needed:

        '3,0'  ~ z^3         O[2,2,2]
        '3,+1' ~ x z^2       O[0,2,2]
        '3,-1' ~ y z^2       O[1,2,2]
        '3,+2' ~ z(x^2-y^2)  O[2,0,0] - O[2,1,1]
        '3,-2' ~ x y z       O[0,1,2]
        '3,+3' ~ x(x^2-3y^2) O[0,0,0] - 3 O[0,1,1]
        '3,-3' ~ y(3x^2-y^2) 3 O[1,0,0] - O[1,1,1]

    Returns
    -------
    (3, 3, 3, d0, d0), fully symmetric and Hermitian in the last two indices.
    It has 7 independent components -- the rank-3 tesseral operators '3,0',
    '3,+-1', '3,+-2', '3,+-3'.

    Returns dict {label: (d0, d0) operator}.
    """

    S = np.asarray(S_low)
    d = S.shape[1]

    M = np.zeros((3, 3, 3, d, d), dtype=complex)
    for a in range(3):
        for b in range(a, 3):
            for c in range(b, 3):

                blk = np.zeros((d, d), dtype=complex)
                for p in permutations((a, b, c)):
                    blk = blk + S[p[0]] @ S[p[1]] @ S[p[2]]
                blk = blk / 6.0

                for p in set(permutations((a, b, c))):
                    M[p] = blk

    V = np.einsum('aabij->bij', M)                       # trace vector
    I3 = np.eye(3)
    octupole_tensor = M - (1 / 5.0) * (np.einsum('ab,cij->abcij', I3, V) +
                                       np.einsum('bc,aij->abcij', I3, V) +
                                       np.einsum('ca,bij->abcij', I3, V))

    x, y, z = 0, 1, 2
    return {'3,0': octupole_tensor[z, z, z],
            '3,+1': octupole_tensor[x, z, z],
            '3,-1': octupole_tensor[y, z, z],
            '3,+2': octupole_tensor[z, x, x] - octupole_tensor[z, y, y],
            '3,-2': octupole_tensor[x, y, z],
            '3,+3': octupole_tensor[x, x, x] - 3 * octupole_tensor[x, y, y],
            '3,-3': 3 * octupole_tensor[y, x, x] - octupole_tensor[y, y, y]}


def determine_phase_rotations(U, S_basis):
    """
    Parameters
    ----------
    U : The rotation that diagonalizes the Sz

    S_basis : Sx, Sy, Sz of a given spin

    Returns
    -------
    U_ph = U*phases; This rotation fixes the Sx, Sy as well as Sz.

    """
    Sc = np.transpose(np.tensordot(U.conj().T, np.tensordot(S_basis, U, axes=(2, 0)), axes=(1, 1)), (1, 0, 2))
    Sp = Sc[0] + 1j * Sc[1]
    ph = np.ones(Sc.shape[1], dtype=complex)
    for i in range(Sc.shape[1] - 1):
        z = Sp[i, i + 1]
        ph[i + 1] = ph[i] * np.conj(z) / abs(z) if abs(z) > 1e-12 else ph[i]
    return U * ph[None, :]


def physical_operator_basis(S_low):
    """
    Parameters
    ----------
    S_low : spin operators

    Description
    -----------
    The raw physical operators -- identity, the three Cartesian pseudospin
    components, and the five tesseral quadrupoles -- built directly from S_low
    with no ladder recursion. Use this basis whenever you need the amplitudes
    attached to specific Cartesian operators.

    Returns
    -------
    ops, ranks (list of ranks of the ops), labels (list of labels of ops),
    scale (list of Frobenius norm of the ops)
    """
    S_low = np.asarray(S_low)
    d = S_low.shape[1]
    ops = [np.eye(d, dtype=complex)]
    nrm = np.sqrt(np.real(np.vdot(ops[0], ops[0])))

    if nrm > 0:
        ops[0] = ops[0] / nrm

    ranks, labels, scale = [0], ['1'], [nrm]

    for a, nm in zip(np.array([1, 2, 0]), np.array(['Sy', 'Sz', 'Sx'])):

        nrm = np.sqrt(np.real(np.vdot(S_low[a], S_low[a])))

        if nrm > 0:
            ops.append(S_low[a] / nrm)
        else:
            ops.append(S_low[a])

        ranks.append(1)
        labels.append(nm)
        scale.append(nrm)

    if d > 2:
        Q = quadrupole_ops(S_low)
        for irrep, Es in quadrupole_pol_basis('spherical').items():

            perm_inds = np.array([2, 4, 0, 3, 1])

            for i, E in enumerate(np.array(Es)[perm_inds]):
                operator = np.einsum('mn,mnij->ij', E, Q)
                nrm = np.sqrt(np.real(np.vdot(operator, operator)))

                if nrm > 0:
                    ops.append(operator / nrm)
                else:
                    ops.append(operator)

                ranks.append(2)
                labels.append(f'Q{i}')
                scale.append(nrm)

    if d > 3:
        perm_inds = np.array([6, 4, 2, 0, 1, 3, 5])

        O_dict = octupole_ops(S_low)

        custom_order = np.array(list(O_dict.keys()))[perm_inds]
        O_reordered_data = {key: O_dict.get(key, None) for key in custom_order}

        for lab, operator in O_reordered_data.items():
            nrm = np.sqrt(np.real(np.vdot(operator, operator)))

            if nrm > 0:
                ops.append(operator / nrm)
            else:
                ops.append(operator)

            ranks.append(3)
            labels.append(f'O{lab[2:]}')
            scale.append(nrm)

    return np.array(ops), np.array(ranks), labels, scale


def stevens_operators(k, spin_ops):
    """
    Parameters
    ----------
    k = rank
    spin_ops = array of spin operators

    Description
    -----------
    dict q -> T^k_q (complex matrix), q = -k..k, built from (S_+)^k by lowering. Standard Racah construction

    Stevens operators : The 2k+1 Hermitian tesseral multipole operators of rank k constructed from dict q-> T^k_q.
    (ops, labels): labels are 'k,+q' (cosine) and 'k,-q' (sine), plus 'k,0'.
      O_k^0    = T^k_0
      O_k^{+q} ~ (T^k_{-q} + (-1)^q T^k_{+q}) / sqrt(2)         (cosine, real)
      O_k^{-q} ~ i (T^k_{-q} - (-1)^q T^k_{+q}) / sqrt(2)       (sine,   real)

    Returns
    -------
    ops, labels

    """

    Sx, Sy, Sz = spin_ops
    Sp = spin_ops[0] + 1j * spin_ops[1]  # S_+ = Sx + i Sy
    Sm = spin_ops[0] - 1j * spin_ops[1]  # S_- = Sx - i Sy
    T = {k: np.linalg.matrix_power(Sp, k).astype(complex)}     # T^k_k ~ (S_+)^k
    for q in range(k, -k, -1):
        T[q - 1] = (-Sm @ T[q] + T[q] @ Sm) / np.sqrt((k + q) * (k - q + 1))

    perm = [k]

    ops, labels = [T[0].copy()], [f'{k},0']
    for q in range(1, k + 1):

        c = ((-1)**q * T[-q] + T[q]) / np.sqrt(2)
        s = 1j * ((-1)**q * T[-q] - T[q]) / np.sqrt(2)

        ops += [c, s]
        labels += [f'{k},+{q}', f'{k},-{q}']
        perm += [k + q, k - q]

    perm_inds = np.argsort(np.array(perm))
    ops = [ops[i] for i in perm_inds]
    labels = [labels[i] for i in perm_inds]

    return ops, labels


def spin_operator_SU2_basis(S_basis, tol=1e-8):
    """
    Parameters
    ----------
    S_basis = exact SU(2) spin operators

    Description
    -----------
    Build Steven operators from exact SU(2) spin operators so the rank d>=2S
    is automatically zero.

    Returns
    -------
    ops, ranks (list of ranks of the ops), labels (list of labels of ops),
    scale (list of Frobenius norm of the ops)

    """
    Sx, Sy, Sz = S_basis
    d = S_basis.shape[1]

    S = (d - 1) / 2
    max_lim = int(2 * S)

    identity = np.eye(d, dtype=complex)

    cand = [(identity, 0, '0,0')]
    for ctr in range(1, max_lim + 1):

        ops, labels = stevens_operators(ctr, S_basis)

        for operator, label in zip(ops, labels):
            cand.append((operator, ctr, label))

    ops, ranks, labels, scale = [], [], [], []

    for operator, k, name in cand:
        nrm = np.sqrt(np.real(np.vdot(operator, operator)))  # Tr(O^dagger O) = Frobenius norm
        if nrm > tol and k < d:                         # keep only operators that exist for this d0

            ops.append(operator / nrm)  # normalise the operator to have unit Frobenius norm
            ranks.append(k)
            labels.append(name)
            scale.append(nrm)

    return np.array(ops), np.array(ranks), labels, scale


def gauge_fix(S_low, tol=1e-10):
    """
    Parameters
    ----------
    S_low : spin operators

    Description
    -----------
    Rotate the multiplet's INTERNAL basis to the standard |S,m> gauge.

    Returns
    -------
    (S_gauged, U) with S_gauged = U^dag S_low U.

    """

    d0 = S_low.shape[1]
    w, U = np.linalg.eigh(S_low[2])
    U = U[:, ::-1]          # m descending

    Sp = np.dot(np.conj(U).T, np.dot((S_low[0] + 1j * S_low[1]), U))

    ph = np.ones(d0, complex)
    for k in range(1, d0):
        z = Sp[k - 1, k]
        ph[k] = ph[k - 1] * np.exp(-1j * np.angle(z)) if abs(z) > tol else ph[k - 1]

    U = U * ph

    return np.einsum('pi,apq,qj->aij', U.conj(), S_low, U), U


def decompose_pseudospin(S_low, tol=1e-8):
    """
    Parameters
    ----------
    S_low : A projected pseudo-spin (low energy levels forming a multiplet)

    Description
    -----------

    Split the projected pseudospin into an exact spin times a
    SPIN-PROJECTION tensor.

        S_low^a  =  sign * sum_b G[a,b] S_exact^b  +  t_a * 1

    with
      G        symmetric positive-definite -- the SPIN-PROJECTION tensor,
               how much of the ideal spin survives projection.
      Rot      the rotation from the lab axes to the principal axes of G
      S_exact  an exact spin-(d0-1)/2 (SU(2) closes to machine precision),
               already rotated into the lab frame: S_exact = Rot . S_ideal
      sign     +1, or -1 when the moment is reversed (normal for a hole)
      t_a      a small trace part, non-zero only because an applied field
               breaks the exact Kramers/time-reversal structure and zero
               otherwise.

    Obtained by least-squares fitting S_low on {1, S_ideal} to get Lambda,
    then polar-decomposing Lambda = G R. Because a rotation preserves SU(2),
    S_exact = R . S_ideal is still an exact spin, and all the anisotropy sits
    in the symmetric factor G.

    Returns
    -------
    dict

    """
    Sin = np.asarray(S_low)
    Sg, Ug = gauge_fix(Sin)

    d0 = Sg.shape[1]
    # only takes Sx, Sy, Sz, not S+/-. These are closed under SU(2) algebra.
    S_id = np.array(spin_ops((d0 - 1) / 2)[:3])

    # forms the basis with the identity and the three spin operators
    basis = [np.eye(d0, dtype=complex)] + list(S_id)

    Gm = np.array([[np.vdot(p, q) for q in basis] for p in basis])  # Gm is diag(d0​, ||S||2, ||S||2, ||S||2)
    # solves the linear system Gm * coef = <p|S[a]> for each a=0,1,2 (x,y,z)
    # to get the coefficients of S_low in the basis of {1, Sx, Sy, Sz}
    coef = np.array([np.linalg.solve(Gm, np.array([np.vdot(p, Sg[a]) for p in basis])) for a in range(3)])
    Lam_c = coef[:, 1:]
    im = float(np.abs(Lam_c.imag).max())
    if im > 1e-8 * max(float(np.abs(Lam_c).max()), 1e-30):
        raise ValueError(f"Lambda has imaginary part {im:.2e} -- is S_low Hermitian?")
    Lam, trace_part = np.real(Lam_c), np.real(coef[:, 0])
    # reconstructs S_low from the coefficients Lam and trace_part
    rec = (np.einsum('ab,bij->aij', Lam, S_id) + np.einsum('a,ij->aij', trace_part, np.eye(d0, dtype=complex)))
    # computes the maximum relative residual between the original S_low and
    # the reconstructed one, normalized by the maximum absolute value of
    # S_low.
    residual = float(np.abs(Sg - rec).max() / max(np.abs(Sg).max(), 1e-30))

    # performs singular value decomposition on Lam to get U, sig (singular
    # values), and Vt. This is used to decompose Lam into a rotation and a
    # symmetric positive-definite matrix.
    U, sig, Vt = np.linalg.svd(Lam)
    Rot = U @ Vt  # computes the rotation matrix Rot from U and Vt
    # determines the sign of the determinant of R, which indicates if the
    # rotation is proper or improper (reflection)
    sign = float(np.sign(np.linalg.det(Rot)))
    if sign < 0:
        Rot = -Rot  # if the determinant is negative, flip the sign of R to make it a proper rotation
    # computes the symmetric positive-definite matrix G from U and the singular values sig
    G = np.dot(U, np.dot(np.diag(sig), U.T))
    # gauged basis, rotates the ideal spin operators S_id by the rotation
    # matrix Rot to get S_exact, which is the exact spin in the lab frame.
    S_ex_g = np.tensordot(Rot, S_id, axes=(1, 0))
    S_exact = np.transpose(
        np.tensordot(
            Ug, np.tensordot(
                S_ex_g, np.conj(Ug).T, axes=(
                    2, 0)), axes=(
                1, 1)), (1, 0, 2))  # back to the original basis
    return dict(G=G, R=Rot, sign=sign, S_exact=S_exact, Lambda=Lam,
                g_principal=sig, trace_part=trace_part, residual=residual)


# ---------------------------------------------------------------------------
#   zero-field splitting -- mandatory for S_pseudo >= 1
# ---------------------------------------------------------------------------
def zfs_tensor(H_low, S_pseudo):
    """
    Parameters
    ----------
    H_low : low energy subspace of the Hamiltonian
    S_pseudo : spin matrices

    Description
    -----------
    D_ab from   P^dag H P - <H>  =  sum_ab D_ab Stilde_a Stilde_b.
    The origin of this term is in spin-orbit coupling and non-degeneracy.

    Returns
    -------
    single-ion anisotropy D matrix

    """
    d0 = H_low.shape[0]
    if d0 < 3:
        return None
    S = (d0 - 1) / 2
    Sv = list(S_pseudo)
    Q = np.array([[0.5 * (Sv[a] @ Sv[b] + Sv[b] @ Sv[a]) - (a == b) * S * (S + 1) / 3 * np.eye(d0)
                 for b in range(3)] for a in range(3)]).reshape(9, d0, d0)
    Hc = H_low - np.trace(H_low) / d0 * np.eye(d0)
    Gm = np.array([[np.vdot(p, q) for q in Q] for p in Q])
    D = np.real(np.linalg.lstsq(Gm, np.array([np.vdot(p, Hc) for p in Q]), rcond=None)[0]).reshape(3, 3)
    D = 0.5 * (D + D.T)

    return D


def odd_polar_blocks(M_low, S_pseudo, tol=1e-10):
    """
    Parameters
    ----------
    M_low     : magnetic moment including lower energy levels
    S_pseudo  : spin-operators exact under SU(2) algebra

    Description
    -----------
    Expand m_low on the T-odd multipoles of every allowed rank (k odd, k <= 2S~).

        blocks[k][a, q] = G^(k)_aq       ops_all[k][q] = O_k^q

    Returns
    -------
    blocks    : coefficients for each Steven operator.
    ops_all   : All the needed odd-rank Steven operators are returned.
    labels, trace_part, reconstructed error, full residual.

    """
    d0 = M_low.shape[1]
    St = (d0 - 1) / 2
    odd = [k for k in range(1, int(2 * St) + 1, 2)]
    identity = np.eye(d0, dtype=complex)

    ops_all, labels = {}, {}
    for k in odd:
        if k == 1:
            ops_all[1] = np.asarray(S_pseudo)          # Unnormalized -> blocks[1] = g
            labels[1] = ['x', 'y', 'z']
        else:
            o, lab = stevens_operators(k, S_pseudo)
            ops_all[k] = np.array([x / np.sqrt(np.real(np.vdot(x, x))) for x in o])
            labels[k] = lab

    def lsq_fit(basis_list):
        B = np.array(basis_list)
        G = np.array([[np.vdot(p, q) for q in B] for p in B])
        W, resid = [], 0.0
        for a in range(3):
            c = np.linalg.solve(G, np.array([np.vdot(p, M_low[a]) for p in B]))
            rec = np.tensordot(c, B, axes=(0, 0))
            resid = max(resid, np.abs(rec - M_low[a]).max() / max(np.abs(M_low[a]).max(), 1e-30))
            W.append(c)
        return np.array(W), resid

    cumulative, resid_ladder = [identity], {}
    for k in odd:                                       # residual after each rank
        cumulative = cumulative + list(ops_all[k])
        W, resid_ladder[k] = lsq_fit(cumulative)      # residual after including rank k
    resid_full = resid_ladder[odd[-1]]

    im = np.abs(W.imag).max()
    if im > 1e-8 * max(np.abs(W).max(), 1e-30):
        raise ValueError(f"Zeeman coefficients have imaginary part {im:.2e} , M_low is not Hermitian")
    W = np.real(W)

    blocks, off = {}, 1
    for k in odd:
        n = len(ops_all[k])
        blocks[k] = W[:, off:off + n]
        off += n
    assert off == W.shape[1], f"block sizes {off} != {W.shape[1]} columns"
    trace_part = W[:, 0]

    rt = reconstruction_error(blocks, ops_all, trace_part, M_low)
    if rt > 1e-10:
        raise ValueError(f"parameter blocks do not rebuild M_low (round trip {rt:.2e})")
    return blocks, ops_all, labels, trace_part, float(rt), resid_full, resid_ladder


def reconstruction_error(blocks, ops_all, trace_part, M_low):
    """
    Parameters
    ----------
    blocks     : Coefficients of ops_all needed to describe the M_low
    ops_all    : All the odd-rank steven operators
    trace_part : The coefficients with identity matrix
    M_low      : Magnetic moment operator

    Description
    -----------
    Resconstructs the M_low = \\sum blocks * ops_all + trace_part * Id

    Returns
    -------
    err        : Error in reconstruction of the M_low

    """
    d0 = M_low.shape[1]
    err = 0.0
    for a in range(3):
        rec = trace_part[a] * np.eye(d0, dtype=complex)
        for k in blocks:
            rec = rec + np.tensordot(blocks[k][a], ops_all[k], axes=(0, 0))
        err = max(err, np.abs(rec - M_low[a]).max() / max(np.abs(M_low[a]).max(), 1e-30))
    return err


def rebuild_moment(blocks, ops_all, trace_part=None, bvec=None):
    """
    Parameters
    ----------
    blocks     :  b_a
    ops_all    :  m_a
    trace_part :  coefficients of identity
    bvec       :  mu_B * (Bx, By, Bz)

    Description
    -----------
    m_a from the reported parameters, or H_Z = sum_a b_a m_a when bvec is
    given (b = mu_B B in eV). It is important use this rather than contracting
    blocks against a single operator list, since rank 1 and ranks >= 3 carry
    different normalisations.

    Returns
    -------
    H_Z

    """
    d0 = next(iter(ops_all.values())).shape[1]
    m = np.zeros((3, d0, d0), dtype=complex)
    for a in range(3):
        if trace_part is not None:
            m[a] = m[a] + trace_part[a] * np.eye(d0, dtype=complex)
        for k in blocks:
            m[a] = m[a] + np.tensordot(blocks[k][a], ops_all[k], axes=(0, 0))
    return m if bvec is None else np.tensordot(np.asarray(bvec), m, axes=(0, 0))


def Multipole_expansion(hmat_i, hmat_n, emat_i_dim, basis_i, trans_ops, l, gamma=0.1,  # noqa: E741
                        ge=2.0, gap_threshold=0.05, gap_tol=1e-9, basis='ideal'):
    """
    Parameters
    ----------
    hmat_i        : Initial Hamiltonian
    hmat_n        : Intermediate Hamiltonian
    emat_i_dim    : dimension of the Valence shell orbitals + Core shell orbitals

    basis_i       : Initial Hamiltonian's Fock basis
    trans_ops     : Transition operators between initial basis to intermediate basis
    l             : The orbital angular momentum of the valence shell;
                    (s, p, d, f) correspond to l = (0, 1, 2, 3)
    gamma         : The inverse-lifetime of the core hole states. Note:
                    currently it is a float but should be an array of floats.
    ge            : Needed to generate the net moment, J = L + ge*S

    gap_threshold : The gap within which the lowest energy levels remain. If a
                    gap keeps a few states from a multiplet, an exception is
                    raised.
    gap_tol       : User-defined tolerance for the difference between the last
                    two energy levels in the multiplet.
    basis         : 'ideal' initiates a Stevens construction of the ladder
                    operators; 'physical' constructs the operators directly.

    Description
    -----------
    Generate the projection operator to the lowest num_spin_states eigenstates
    of the initial Hamiltonian. It projects the RIXS cross-section to the
    low-energy subspace of the initial Hamiltonian.
    gap_threshold  # eV: >> ext_B splitting, << first dd/exchange gap

    Returns
    -------
    dict

    """

    eval_i, evec_i = np.linalg.eigh(hmat_i)
    E0 = eval_i[0]

    rel = eval_i - E0

    # The low energy subspace #

    d0_raw = max(int(np.searchsorted(rel, gap_threshold)), 2)
    d0 = d0_raw
    while d0 < len(rel) and rel[d0] - rel[d0 - 1] < gap_tol:
        d0 += 1

    if d0 != d0_raw:
        raise Exception("The gap threshold is chosen such that it slices a multiplet, increase the threshold to " +
                        str(rel[d0] if d0 < len(rel) else np.inf))

    # P is the projection operator to the lowest num_spin_states eigenstates of the initial Hamiltonian.
    P = evec_i[:, :d0]

    # Determine the pseudo-spin

    L_d = edrixs.get_orb_momentum(l, ispin=True)
    S_d = edrixs.get_spin_momentum(l)
    M_d = L_d + ge * S_d
    # Embed the 10x10 d-spin operators into the FULL single-particle space (valence d + core p = 16 orbitals).
    M_spin_single_particle = np.zeros((3, emat_i_dim, emat_i_dim), dtype=complex)
    M_spin_single_particle[:, :M_d.shape[1], :M_d.shape[2]] = M_d

    M_spin = np.stack([edrixs.build_op(component, None, lb=basis_i, rb=basis_i, backend="dense")
                      for component in M_spin_single_particle])

    M_low = np.transpose(np.tensordot(np.tensordot(np.conj(P).T, M_spin, axes=(1, 1)), P, axes=(2, 0)), (1, 0, 2))

    dec = decompose_pseudospin(M_low)
    g_mat = dec['sign'] * dec['G']
    GG = g_mat @ g_mat.T

    S_pseudo = dec['S_exact']             # USE THIS pseudospin everywhere
    residual = dec['residual']

    H_low = np.dot(np.conj(P).T, np.dot(hmat_i, P))
    D = zfs_tensor(H_low, S_pseudo)  # determines the single-ion anisotropy or zero field splitting term

    blocks, ops_all, labels_G, trace_part, rt, resid_full, resid_ladder = odd_polar_blocks(M_low, S_pseudo)

    dg = np.abs(blocks[1] - g_mat).max()
    if dg > 1e-8:
        raise ValueError(f"blocks[1] is not the g-matrix (max diff {dg:.2e})")

    # ---- RIXS operator basis ----------------------------------------------
    if basis == 'ideal':
        B_ops, ranks, labels, scale = spin_operator_SU2_basis(S_pseudo)
    else:
        B_ops, ranks, labels, scale = physical_operator_basis(S_pseudo)

    # diagonalize the intermediate Hamiltonian to get eigenvalues and
    # eigenvectors, assumes intermediate Hamiltonian is Hermitian.
    e_val_n, e_vec_n = np.linalg.eigh(hmat_n)
    # transform the transition operators to the eigenbasis of the intermediate Hamiltonian.
    trans_ops_mod = np.tensordot(np.conj(e_vec_n).T, trans_ops, axes=(1, 1))
    trans_ops_mod = np.transpose(trans_ops_mod, (1, 0, 2))

    return dict(
        # manifold
        d0=d0, Stilde=(d0 - 1) / 2, P=P, E0=E0, rel=rel, H_low=H_low,
        # pseudospin and Zeeman parameters
        S_basis=S_pseudo,          # unnormalised, ONE convention
        g_matrix=g_mat, GG=GG, residual=residual, roundtrip=rt,
        resid_full=resid_full, resid_ladder=resid_ladder,
        blocks=blocks, ops_all=ops_all, block_labels=labels_G, trace_part=trace_part, D=D, M_low=M_low,
        # RIXS
        B_ops=B_ops, ranks=ranks, labels=labels, scale=scale,
        trans_ops_mod=trans_ops_mod, e_val_n=e_val_n,
    )


def resonant_omegas(eval_n, E0, trans_ops_eig, P, gamma, rel_thresh=1e-2, cluster_tol=None, eps_in=None):
    """
    Parameters
    ----------
    eval_n        :    Intermediate Hamiltonian's eigenvalues.
    E0            :    Ground state energy of the initial Hamiltonian.
    trans_ops_eig :    The transition operators from initial to intermediate
                       Hamiltonian in the initial and intermediate Fock basis.
    P             :    Projection operator from the full basis (valence shell) to lower energy basis.
    gamma         :    The inverse-lifetime of the core hole states. Currently
                       a float, but it should be an array of floats.
    rel_thresh    :    Threshold of the XAS intensity of the resonant peaks. Default is 1%.
    cluster_tol   :    Spread of peaks within this range is grouped in a single peak.
    eps_in        :    Polarization of the incoming x-rays.

    Description
    -----------
    Poles of the KH resolvent along omega_in sit at eval_n - E0. Keep the most
    spectrally intense (bright) ones (absorption weight > rel_thresh*max) and
    merge poles closer than cluster_tol (default gamma) into
    intensity-weighted centroids. By default we keep peaks with 1% of the
    maximum absorption weight and merge peaks closer than the core-hole
    lifetime broadening gamma.

    Returns
    -------
    omega_res : Array of resonant omega values
    w_res     : Intensity of the RIXS in omega values

    """

    # <c| D_m |g>, (3, N_n, d0) Alternatively Aabs = np.einsum("mci,ig->mcg", trans_ops_eig, P)
    Aabs = np.tensordot(trans_ops_eig, P, axes=(2, 0))
    if eps_in is None:
        xas_w = np.sum(np.abs(Aabs)**2, axis=(0, 2))        # isotropic (all polarizations)
    else:
        xas_w = np.sum(np.abs(np.tensordot(eps_in, Aabs, axes=(0, 0)))**2, axis=1)

    res = eval_n - E0

    if xas_w.max() <= 0:
        return np.array([]), np.array([])

    keep = xas_w > rel_thresh * xas_w.max()
    res_k, w_k = res[keep], xas_w[keep]
    order = np.argsort(res_k)
    res_k, w_k = res_k[order], w_k[order]

    if cluster_tol is None:
        cluster_tol = gamma

    split = np.where(np.diff(res_k) > cluster_tol)[0] + 1
    groups = np.split(np.arange(len(res_k)), split)

    omega_res = np.array([np.sum(w_k[g] * res_k[g]) / np.sum(w_k[g]) for g in groups])
    w_res = np.array([np.sum(w_k[g]) for g in groups])

    return omega_res, w_res


def generate_rixs_coefficients(omega_res, e_val_n, E0, trans_ops_mod, P, gamma, B_ops, verbose=False):
    """
    Parameters
    ----------
    omega_res     :     List of omegas for which RIXS operators are computed.
    e_val_n       :     Intermediate Hamiltonian's eigenvalues.
    E0            :     Ground state energy of the initial Hamiltonian.
    trans_ops_mod :     Transition operators from initial to intermediate
                        Hamiltonian eigenstates (U_n^\\dag T_q).
    P             :     Projection operator from the full basis (valence shell) to lower energy basis.
    gamma         :     The inverse-lifetime of the core hole states. Currently
                        a float, but it should be an array of floats.
    B_ops         :     The basis of stevens operators.
    verbose       :     'True' prints all the resonant omega values. Default is 'False'.

    Raises
    ------
    an exception is the Basis reconstruction is incomplete then R is not reconstructed.

    Returns
    -------
    C             :     Coefficients for each operator in the basis B_ops such
                        that R[w] = \\sum C[w]*B_ops, where w is a resonant
                        omega value.
    R_list        :     List of R[w] operators at each resonant energy values.

    """

    n_ops = B_ops.shape[0]
    C = np.zeros((len(omega_res), 3, 3, n_ops), dtype=complex)

    # Gram matrix of the basis (exact projection needs it; near-identity here).
    Gram = np.tensordot(B_ops.conj(), B_ops, axes=([1, 2], [1, 2]))
    Gram_inv = np.linalg.inv(Gram)

    R_list = []
    for w, omega_in in enumerate(omega_res):
        if verbose:
            print(f"Resonance no. {w + 1}: {omega_in:.2f} eV")                 # all the edges -> RIXS map
        z = E0 + omega_in + 1j * gamma / 2
        # evaluate the denominator of the Kramers-Heisenberg formula, which is the
        # Green's function of the intermediate state Hamiltonian.
        Green = np.diag(1.0 / (e_val_n - z))
        R_valence = np.tensordot(
            np.tensordot(
                np.conj(
                    np.transpose(
                        trans_ops_mod, (0, 2, 1))), Green, axes=(
                    2, 0)), trans_ops_mod, axes=(
                        2, 1))
        R_valence = np.transpose(R_valence, (0, 2, 1, 3))

        R = np.tensordot(np.tensordot(np.conj(P).T, R_valence, axes=(1, 2)), P, axes=(3, 0))
        R = np.transpose(R, (1, 2, 0, 3))
        R_list.append(R)

        M = np.transpose(np.tensordot(np.conj(B_ops), R, axes=([1, 2], [2, 3])), (1, 2, 0))
        C[w, :, :, :] = np.tensordot(M, Gram_inv.T, axes=(2, 0))  # (m,n,a) = sum_b M[m,n,b] * Gram_inv[b,a]

        err = np.linalg.norm(np.tensordot(C[w], B_ops, axes=(2, 0)) - R)
        if err > 1e-8:
            raise ValueError(f"basis incomplete: reconstruction error {err:.2e}")

    return C, R_list


if __name__ == "__main__":

    # Ni2+ ion in La2NiO4 compound from Fabbris et al., PRL 118, 156402
    # (2017). Scaling percentage taken from the paper and actual values are
    # obtained from get_atom_data.py script in the edrixs package.

    # Number of occupancy of 3d shell
    valence_noccu = 8

    shell_level_energies = (0, 0)

    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=valence_noccu, edge='L23')
    name_i, slat_i = [list(i) for i in zip(*res['slater_i'])]
    name_n, slat_n = [list(i) for i in zip(*res['slater_n'])]

    # Slater integrals for initial Hamiltonian without core-hole
    slater_i = edrixs.rescale(slat_i, ([1, 2], [0.65, 0.65]))
    slater_i[0] = edrixs.get_F0('d', slater_i[1], slater_i[2])    # F0_dd

    # Slater integrals for intermediate Hamiltonian with core-hole
    slater_n = edrixs.rescale(slat_n, ([1, 2, 4, 5, 6], [0.65, 0.65, 0.65, 0.85, 0.85]))
    slater_n[0] = edrixs.get_F0('d', slater_n[1], slater_n[2])     # F0_dd
    slater_n[3] = edrixs.get_F0('dp', slater_n[5], slater_n[6])    # F0_dp

    slater = (slater_i, slater_n)

    # Spin-orbit coupling strengths
    zeta_d_i = res['v_soc_i'][0]  # valence 3d electron without core-hole
    zeta_d_n = res['v_soc_n'][0]  # valence 3d electron with core-hole
    # E_{L2} - E_{L3} = 1.5 * zeta_p
    zeta_p_n = (res['edge_ene'][0] - res['edge_ene'][1]) / 1.5  # core 2p electron
    valence_spin_orbital_coupling = (zeta_d_i, zeta_d_n)
    core_spin_orbital_coupling = zeta_p_n

    v_cfmat = edrixs.cf_tetragonal_d(1.6, 0.75, 0.1)
    name = "Ni2_ion_D4h_symmetry_La2NiO4"

    # inverse core-hole lifetime broadening (in eV) for the intermediate Hamiltonian.
    gamma = 0.275

    # theta_in + theta_out equals the angle between k_in and k_out

    two_theta = (150 / 180) * np.pi                          # total scattering angle
    theta_in = (10 / 180) * np.pi                           # grazing incidence angle (pi/2 = normal incidence)
    theta_out = two_theta - theta_in                     # edrixs: theta_in + theta_out = 2*theta
    phi = 0.0

    incident_pol = 'pi'  # 'sigma' or 'pi'
    scattered_pol = 'pi'   # 'sigma' or 'pi'

    out = edrixs.model_1v1c(shell_name=('d', 'p'), shell_level=shell_level_energies,
                            v_soc=valence_spin_orbital_coupling,
                            c_soc=core_spin_orbital_coupling,
                            v_noccu=valence_noccu, slater=(slater_i, slater_n),
                            ext_B=[0, 0, 0.0], on_which='both',
                            v_cfmat=v_cfmat)

    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat = out

    hmat_i, hmat_n, trans_ops = edrixs.get_ops(
        emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat, backend="dense")

    trans_ops = np.asarray(trans_ops)

    emat_i_dim = emat_i.shape[0]

    orbital_angular_momentum = 2  # d shell (s:0, p:1, d:2, f:3)

    G_E = 2.0
    gap_threshold = 0.05
    gap_tol = 1e-9
    basis = 'ideal'

    output_dict = Multipole_expansion(
        hmat_i,
        hmat_n,
        emat_i_dim,
        basis_i,
        trans_ops,
        orbital_angular_momentum,
        gamma,
        G_E,
        gap_threshold,
        gap_tol,
        basis)

    num_spin_states = output_dict["d0"]
    B_ops = output_dict["B_ops"]
    ranks = output_dict["ranks"]
    labels = output_dict["labels"]
    S_basis = output_dict["S_basis"]
    P = output_dict["P"]
    D = output_dict["D"]
    residual = output_dict["residual"]
    scale = output_dict["scale"]
    G_mat = output_dict["g_matrix"]
    trans_ops_mod = output_dict["trans_ops_mod"]
    e_val_n = output_dict["e_val_n"]
    E0 = output_dict["E0"]
    resid_full = output_dict["resid_full"]
    resid_ladder = output_dict["resid_ladder"]

    print("The residual in fitting the rank-1 operators with full moment is " + str(resid_ladder[1]))
    if resid_ladder[1] > 1e-6:
        print("The full residual is " + str(resid_full) + " after including odd ranked operators upto rank " +
              str(S_basis.shape[1] if S_basis.shape[1] % 2 == 1 else S_basis.shape[1] - 1))

    # Polarization vectors for RIXS

    CHAN = {'pi-pi': (0.0, 0.0),
            'sigma-pi': (np.pi / 2, 0.0),
            'pi-sigma': (0.0, np.pi / 2),
            'sigma-sigma': (np.pi / 2, np.pi / 2)}

    if incident_pol == 'sigma' and scattered_pol == 'sigma':
        eps_in, eps_out = edrixs.dipole_polvec_rixs(theta_in, theta_out, phi, *CHAN['sigma-sigma'])
        str_pol = "sigma_in_sigma_out"
    if incident_pol == 'sigma' and scattered_pol == 'pi':
        eps_in, eps_out = edrixs.dipole_polvec_rixs(theta_in, theta_out, phi, *CHAN['sigma-pi'])
        str_pol = "sigma_in_pi_out"
    if incident_pol == 'pi' and scattered_pol == 'sigma':
        eps_in, eps_out = edrixs.dipole_polvec_rixs(theta_in, theta_out, phi, *CHAN['pi-sigma'])
        str_pol = "pi_in_sigma_out"
    if incident_pol == 'pi' and scattered_pol == 'pi':
        eps_in, eps_out = edrixs.dipole_polvec_rixs(theta_in, theta_out, phi, *CHAN['pi-pi'])
        str_pol = "pi_in_pi_out"

    # run RIXS only at the bright resonances of isotropic XAS.
    omega_res, w_res = resonant_omegas(e_val_n, E0, trans_ops_mod, P, gamma, eps_in=eps_in)

    # generate the RIXS coefficients for each resonance and each polarization pair.
    C, R_list = generate_rixs_coefficients(omega_res, e_val_n, E0, trans_ops_mod, P, gamma, B_ops, verbose=True)

    comm = S_basis[0] @ S_basis[1] - S_basis[1] @ S_basis[0]  # the commutator of Sx and Sy
    assert np.abs(comm - 1j * S_basis[2]).max() < 1e-10, "S_basis does not close into SU(2)"

    # Convert from exact SU(2) operators to the usual spin operators

    sunny = True
    str_sunny = "sunny_" if sunny else "unit_Frobenius_norm_"

    if sunny:
        c_kq = np.ones(len(labels))
        for i, lab in enumerate(labels):
            if ',' not in lab:                     # 'Sx', 'Q0', '1' -- physical basis
                continue
            k, q = (int(t) for t in lab.split(','))
            if k == 0:                             # the identity
                continue
            c_kq[i] = stevens_to_sunny(k, q)
        s_scale = np.array(scale) / c_kq             # = ||O_Sunny||_F   (1 for rank 0)
        B_ops = s_scale[:, None, None] * B_ops     # no [i]
        C = C / s_scale[None, None, None, :]

    _, vec = np.linalg.eigh(S_basis[2])
    U_rot_cart = determine_phase_rotations(vec[:, ::-1], S_basis)

    S_cart = np.asarray([np.dot(U_rot_cart.conj().T, np.dot(S_basis[ctr], U_rot_cart)) for ctr in range(3)])

    R_cart_list = np.zeros_like(np.asarray(R_list), dtype=complex)
    B_ops_cart = np.zeros_like(B_ops, dtype=complex)

    for i in range(len(R_list)):
        R_cart_list[i] = np.tensordot(
            U_rot_cart.conj().T, np.tensordot(
                R_list[i], U_rot_cart, axes=(
                    3, 0)), axes=(
                1, 2))
        R_cart_list[i] = np.transpose(R_cart_list[i], (1, 2, 0, 3))

    for i in range(len(B_ops)):
        B_ops_cart[i] = np.dot(np.conj(U_rot_cart).T, np.dot(B_ops[i], U_rot_cart))

    hdf5_file_path = f'spin_matrices_{name}_polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('spin_matrices', data=S_cart)
    hdf5_file_path = f'R_list_{name}_{str_sunny}polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('R_list', data=np.asarray(R_cart_list))
    hdf5_file_path = f'Steven_ops_basis_{name}_{str_sunny}polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('Steven_op_basis', data=B_ops_cart)
    hdf5_file_path = f'Coefficients_Steven_ops_basis_{name}_{str_sunny}polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('Coeff_Steven_op_basis', data=C)

    hdf5_file_path = f'gmatrix_{name}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('gmatrix', data=G_mat)
    if num_spin_states >= 3:
        hdf5_file_path = f'single_ion_anisotropy_D_{name}.h5'
        with h5py.File(hdf5_file_path, 'w') as hdf5_file:
            hdf5_file.create_dataset('D', data=D)
    hdf5_file_path = f'resonant_omegas_{name}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('resonant_omegas', data=omega_res)

    hdf5_file_path = f'polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('pol_vecs_epsin_epsout', data=np.array([eps_in, eps_out]))

    # picking the first resonance for the RIXS cross-section, which is the lowest energy resonance.
    w = 0

    A = np.tensordot(eps_out.conj(), np.tensordot(C[w], eps_in, axes=(1, 0)), axes=(0, 0))

    hdf5_file_path = f'A_Coefficients_Steven_ops_basis_{name}_{str_sunny}polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('Coeff_Steven_op_basis', data=A)

    c0 = A[ranks == 0][0]                      # identity amplitude (charge / elastic, non-spin-flip)
    print("c0 (identity / elastic amplitude):", np.round(c0, 6))

    A_spin = A[ranks == 1]                     # (Ax, Ay, Az) -- the single-magnon amplitudes

    # Single-ion magnetic scattering tensor F_ab = A_a A_b^*.
    # The single-magnon RIXS intensity on a lattice is  I(q,w) ~ sum_ab F_ab S_ab(q,w).
    F = np.outer(A_spin, A_spin.conj())

    df = pd.DataFrame(F, columns=['Sy', 'Sz', 'Sx'], index=['Sy', 'Sz', 'Sx'])
    print(df)

    if num_spin_states > 2:
        A_quadrupole = A[ranks == 2]           # (Q1, Q2, Q3, Q4, Q5) -- the two-magnon amplitudes
        F_quad = np.outer(A_quadrupole, A_quadrupole.conj())
        df = pd.DataFrame(
            F_quad,
            columns=[
                '{Sx, Sy}',
                '{Sz, Sy}',
                '3Sz2 - S(S+1)',
                '{Sz, Sx}',
                'Sx2 - Sy2'],
            index=[
                '{Sx, Sy}',
                '{Sz, Sy}',
                '3Sz2 - S(S+1)',
                '{Sz, Sx}',
                'Sx2 - Sy2'])
        print(df)

    if num_spin_states > 3:
        A_octupole = A[ranks == 3]             # (O1, O2, O3, O4, O5, O6, O7) -- the three-magnon amplitudes
        F_oct = np.outer(A_octupole, A_octupole.conj())
        df = pd.DataFrame(
            F_oct,
            columns=[
                'y(3x2-y2)',
                'xyz',
                'y(5z2 - r2)',
                'z(5z2 - 3r2)',
                'x(5z2 - r2)',
                'z(x2-y2)',
                'x(x2 - 3y2)'],
            index=[
                'y(3x2-y2)',
                'xyz',
                'y(5z2 - r2)',
                'z(5z2 - 3r2)',
                'x(5z2 - r2)',
                'z(x2-y2)',
                'x(x2 - 3y2)'])
        print(df)

    F = np.outer(A, A.conj())

    hdf5_file_path = f'F_{name}_{str_sunny}polarization_{str_pol}.h5'
    with h5py.File(hdf5_file_path, 'w') as hdf5_file:
        hdf5_file.create_dataset('F_mat', data=F)

    df = pd.DataFrame(F, columns=labels, index=labels)  # works for any S
    print(df.round(7))

    print(np.dot(G_mat.T, G_mat))

    print(D)
