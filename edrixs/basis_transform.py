__all__ = ['cb_op', 'cb_op2', 'tmat_c2r', 'tmat_r2c', 'tmat_r2cub_f',
           'tmat_cub2r_f', 'tmat_c2j', 'tmat_r2trig', 'tmat_trig2r',
           'transform_utensor', 'fourier_hr2hk']

import numpy as np


def cb_op(oper_O, TL, TR=None):
    """
    Change the basis of an operator :math:`\\hat{O}`.

    .. math::

        O^{\\prime} = (T_L)^{\\dagger} O (T_R),

    Parameters
    ----------
    oper_O: array-like
        At least 2-dimension, the last 2-dimension is the dimension of the
        matrix form of operator :math:`\\hat{O}` in basis :math:`A`.
        For example:

        - oper_O.shape = (3, 10, 10), means 3 operators with dimension :math:`10 \\times 10`

        - oper_O.shape = (2, 3, 10, 10), means :math:`2 \\times 3=6` operators with
          dimension :math:`10 \\times 10`
    TL: 2d array
        The unitary transformation matrix from basis :math:`A` to
        basis :math:`B`,  namely,

        :math:`TL_{ij} = <\\psi^{A}_{i}|\\phi^{B}_{j}>`.

    TR: 2d array
        The unitary transformation matrix from basis :math:`A` to
        basis :math:`B`,  namely,

        :math:`TR_{ij} = <\\psi^{A}_{i}|\\phi^{B}_{j}>`.

        if TR = None, TR = TL

    Returns
    -------
    res: same shape as oper_O
        The matrices form of operators :math:`\\hat{O}` in new basis.
    """
    oper_O = np.array(oper_O, order='C')
    dim = oper_O.shape
    if TR is None:
        TR = TL
    if len(dim) < 2:
        raise Exception("Dimension of oper_O should be at least 2")
    elif len(dim) == 2:
        res = np.dot(np.dot(np.conj(np.transpose(TL)), oper_O), TR)
    else:
        tot = np.prod(dim[0:-2])
        tmp_oper = oper_O.reshape((tot, dim[-2], dim[-1]))
        for i in range(tot):
            tmp_oper[i] = np.dot(np.dot(np.conj(np.transpose(TL)), tmp_oper[i]), TR)
        res = tmp_oper.reshape(dim)

    return res


def cb_op2(oper_O, TL, TR):
    """
    Change the basis of an operator :math:`\\hat{O}`.

    .. math::

        O^{\\prime} = (T_L)^{\\dagger} O (T_R),


    Parameters
    ----------
    oper_O: array-like
        At least 2-dimension, the last 2-dimension is the dimension of the
        matrix form of operator :math:`\\hat{O}` in basis :math:`A`.
        For example:

        - oper_O.shape = (3, 10, 10), means 3 operatos with dimension :math:`10 \\times 10`

        - oper_O.shape = (2, 3, 10, 10), means :math:`2 \\times 3=6` operators with
          dimension :math:`10 \\times 10`
    TL: 2d array
        The unitary transformation matrix applied on the left.
    TR: 2d array
        The unitary transformation matrix applied on the right.

    Returns
    -------
    res: same shape as oper_O
        The matrices form of operators :math:`\\hat{O}` in new basis.
    """
    oper_O = np.array(oper_O, order='C')
    dim = oper_O.shape
    if len(dim) < 2:
        raise Exception("Dimension of oper_O should be at least 2")
    elif len(dim) == 2:
        res = np.dot(np.dot(np.conj(np.transpose(TL)), oper_O), TR)
    else:
        tot = np.prod(dim[0:-2])
        tmp_oper = oper_O.reshape((tot, dim[-2], dim[-1]))
        for i in range(tot):
            tmp_oper[i] = np.dot(np.dot(np.conj(np.transpose(TL)), tmp_oper[i]), TR)
        res = tmp_oper.reshape(dim)

    return res


def tmat_c2r(case, ispin=False):
    """
    Get the unitary transformation matrix from the basis of complex
    spherical harmonics to real spherical harmonics.

    Parameters
    ----------
    case: string
        Label for different systems.

        - 'p': for :math:`p`-shell
        - 't2g': for :math:`t_{2g}`-shell
        - 'd': for :math:`d`-shell
        - 'f': for :math:`f`-shell
    ispin: logical
        Whether to include spin degree of freedom or not (default: False).

    Returns
    -------
    t_c2r: 2d complex array
        The transformation matrix.
    """

    sqrt2 = np.sqrt(2.0)
    ci = np.complex128(0.0 + 1.0j)
    cone = np.complex128(1.0 + 0.0j)

    # p orbitals: px, py, pz
    if case.strip() == 'p':
        norbs = 3
        t_c2r = np.zeros((norbs, norbs), dtype=np.complex128)
        # px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0, 0] = cone / sqrt2
        t_c2r[2, 0] = -cone / sqrt2
        # py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0, 1] = ci / sqrt2
        t_c2r[2, 1] = ci / sqrt2
        # pz=|1,0>
        t_c2r[1, 2] = cone

    # t2g orbitals in the t2g subspace, here, we use the so-called
    # T-P equivalence, t2g orbitals behave like the effective orbital
    # angular momentum leff=1
    # dzx ~ py,  dzy ~ px, dxy ~ pz
    elif case.strip() == 't2g':
        norbs = 3
        t_c2r = np.zeros((norbs, norbs), dtype=np.complex128)
        # dzx --> py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0, 0] = ci / sqrt2
        t_c2r[2, 0] = ci / sqrt2
        # dzy --> px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0, 1] = cone / sqrt2
        t_c2r[2, 1] = -cone / sqrt2
        # dxy --> pz=|1,0>
        t_c2r[1, 2] = cone

    # d orbitals: dz2, dzx, dzy, dx2-y2, dxy
    elif case.strip() == 'd':
        norbs = 5
        t_c2r = np.zeros((norbs, norbs), dtype=np.complex128)
        # dz2=|2,0>
        t_c2r[2, 0] = cone
        # dzx=1/sqrt(2)( |2,-1> - |2,1> )
        t_c2r[1, 1] = cone / sqrt2
        t_c2r[3, 1] = -cone / sqrt2
        # dzy=i/sqrt(2)( |2,-1> + |2,1> )
        t_c2r[1, 2] = ci / sqrt2
        t_c2r[3, 2] = ci / sqrt2
        # dx2-y2=1/sqrt(2)( |2,-2> + |2,2> )
        t_c2r[0, 3] = cone / sqrt2
        t_c2r[4, 3] = cone / sqrt2
        # dxy=i/sqrt(2)( |2,-2> - |2,2> )
        t_c2r[0, 4] = ci / sqrt2
        t_c2r[4, 4] = -ci / sqrt2

    # f orbitals, please NOTE that this real form of the f orbitals is not the
    # basis of the representation of the cubic point group, please call the
    # function ``tmat_r2cub" to get the transformation matrix from this basis
    # to the cubic basis that is the representation of the cubic point group.
    elif case.strip() == 'f':
        norbs = 7
        t_c2r = np.zeros((norbs, norbs), dtype=np.complex128)
        # fz3 = |3,0>
        t_c2r[3, 0] = cone
        # fxz2 = 1/sqrt(2)( |3,-1> - |3,1> )
        t_c2r[2, 1] = cone / sqrt2
        t_c2r[4, 1] = -cone / sqrt2
        # fyz2 = i/sqrt(2)( |3,-1> + |3,1> )
        t_c2r[2, 2] = ci / sqrt2
        t_c2r[4, 2] = ci / sqrt2
        # fz(x2-y2) = 1/sqrt(2)( |3,-2> + |3,2> )
        t_c2r[1, 3] = cone / sqrt2
        t_c2r[5, 3] = cone / sqrt2
        # fxyz = i/sqrt(2)( |3,-2> - |3,2> )
        t_c2r[1, 4] = ci / sqrt2
        t_c2r[5, 4] = -ci / sqrt2
        # fx(x2-3y2) = 1/sqrt(2) ( |3,-3> - |3,3> )
        t_c2r[0, 5] = cone / sqrt2
        t_c2r[6, 5] = -cone / sqrt2
        # fy(3x2-y2) = i/sqrt(2) ( |3,-3> + |3,3> )
        t_c2r[0, 6] = ci / sqrt2
        t_c2r[6, 6] = ci / sqrt2
    else:
        raise Exception("error in tmat_c2r: Do NOT support tmat_c2r for this case: ", case)

    # the spin order is: up dn up dn ... up dn
    if ispin:
        ntot_orbs = 2 * norbs
        t_c2r_spin = np.zeros((ntot_orbs, ntot_orbs), dtype=np.complex128)
        # spin up
        t_c2r_spin[0:ntot_orbs:2, 0:ntot_orbs:2] = t_c2r
        # spin dn
        t_c2r_spin[1:ntot_orbs:2, 1:ntot_orbs:2] = t_c2r
        return t_c2r_spin
    else:
        return t_c2r


def tmat_r2c(case, ispin=False):
    """
    Get the unitary transformation matrix from the basis of real
    spherical harmonics to complex spherical harmonics.

    Parameters
    ----------
    case: string
        Label for different systems.

        - 'p': for :math:`p`-shell
        - 't2g': for :math:`t_{2g}`-shell
        - 'd': for :math:`d`-shell
        - 'f': for :math:`f`-shell
    ispin: logical
        Whether to include spin degree of freedom or not (default: False).

    Returns
    -------
    t_r2c: 2d complex array
        The transformation matrix.
    """

    t_r2c = np.conj(np.transpose(tmat_c2r(case, ispin)))
    return t_r2c


def tmat_r2cub_f(ispin=False):
    """
    Get the transformation matrix from real spherical harmonics to the
    cubic spherical harmonics that is the representation of the cubic
    point group, only for :math:`f` system.

    Parameters
    ----------
    ispin: logical
        Whether including spin degree of freedom or not (default: False).

    Returns
    -------
    t_r2cub: 2d complex array
        The transformation matrix.
    """

    a = np.sqrt(10.0) / 4.0 + 0.0j
    b = np.sqrt(6.0) / 4.0 + 0.0j
    c = 1.0 + 0.0j

    norbs = 7
    t_r2cub = np.zeros((norbs, norbs), dtype=np.complex128)

    # T1u
    # fx3 = -sqrt(6)/4 fxz2 + sqrt(10)/4 fx(x2-3y2)
    t_r2cub[1, 0] = -b
    t_r2cub[5, 0] = a
    # fy3 = -sqrt(6)/4 fyz2 - sqrt(10)/4 fy(3x2-y2)
    t_r2cub[2, 1] = -b
    t_r2cub[6, 1] = -a
    # fz3 = fz3
    t_r2cub[0, 2] = c

    # T2u
    # fx(y2-z2) = -sqrt(10)/4 fxz2 - sqrt(6)/4 fx(x2-3y2)
    t_r2cub[1, 3] = -a
    t_r2cub[5, 3] = -b
    # fy(z2-x2) = sqrt(10)/4 fyz2 - sqrt(6)/4 fy(3x2-y2)
    t_r2cub[2, 4] = a
    t_r2cub[6, 4] = -b
    # fz(x2-y2) = fz(x2-y2)
    t_r2cub[3, 5] = c

    # A2u
    # fxyz = fxyz
    t_r2cub[4, 6] = c

    if ispin:
        ntot_orbs = 2 * norbs
        t_r2cub_spin = np.zeros((ntot_orbs, ntot_orbs), dtype=np.complex128)
        # spin up
        t_r2cub_spin[0:ntot_orbs:2, 0:ntot_orbs:2] = t_r2cub
        # spin dn
        t_r2cub_spin[1:ntot_orbs:2, 1:ntot_orbs:2] = t_r2cub
        return t_r2cub_spin
    else:
        return t_r2cub


def tmat_cub2r_f(ispin=False):
    """
    Get the transformation matrix from the cubic spherical harmonics to
    real spherical harmonics, only for :math:`f` system.

    Parameters
    ----------
    ispin: logical
        Whether including spin degree of freedom or not (default: False).

    Returns
    -------
    t_cub2r: 2d complex array
        The transformation matrix.
    """

    t_cub2r = np.conj(np.transpose(tmat_r2cub_f(ispin)))
    return t_cub2r


def tmat_c2j(orb_l):
    """
    Get the transformation matrix from the complex spherical harmonics to
    the :math:`|j^2,j_z>` basis in which the spin-oribt coupling Hamiltonian
    is diagonal. The orbital order is:

    :math:`|j=l-1/2, -j>, |j=l-1/2, -j+1>, ... |j=l-1/2, +j>,`

    :math:`|j=l+1/2, -j>, |j=l+1/2, -j+1>, ..., |j=l+1/2, +j>`.

    Parameters
    ----------
    orb_l: int
        Quantum number of orbital angular momentum.

    Returns
    -------
    t_c2j: 2d complex array
        The transformation matrix.
    """

    if orb_l == 1:
        t_c2j = np.zeros((6, 6), dtype=np.complex128)
        t_c2j[0, 0] = -np.sqrt(2.0 / 3.0)
        t_c2j[3, 0] = np.sqrt(1.0 / 3.0)
        t_c2j[2, 1] = -np.sqrt(1.0 / 3.0)
        t_c2j[5, 1] = np.sqrt(2.0 / 3.0)
        t_c2j[1, 2] = 1.0
        t_c2j[0, 3] = np.sqrt(1.0 / 3.0)
        t_c2j[3, 3] = np.sqrt(2.0 / 3.0)
        t_c2j[2, 4] = np.sqrt(2.0 / 3.0)
        t_c2j[5, 4] = np.sqrt(1.0 / 3.0)
        t_c2j[4, 5] = 1.0

        return t_c2j

    elif orb_l == 2:
        t_c2j = np.zeros((10, 10), dtype=np.complex128)
        t_c2j[0, 0] = -np.sqrt(4.0 / 5.0)
        t_c2j[3, 0] = np.sqrt(1.0 / 5.0)
        t_c2j[2, 1] = -np.sqrt(3.0 / 5.0)
        t_c2j[5, 1] = np.sqrt(2.0 / 5.0)
        t_c2j[4, 2] = -np.sqrt(2.0 / 5.0)
        t_c2j[7, 2] = np.sqrt(3.0 / 5.0)
        t_c2j[6, 3] = -np.sqrt(1.0 / 5.0)
        t_c2j[9, 3] = np.sqrt(4.0 / 5.0)
        t_c2j[1, 4] = 1.0
        t_c2j[0, 5] = np.sqrt(1.0 / 5.0)
        t_c2j[3, 5] = np.sqrt(4.0 / 5.0)
        t_c2j[2, 6] = np.sqrt(2.0 / 5.0)
        t_c2j[5, 6] = np.sqrt(3.0 / 5.0)
        t_c2j[4, 7] = np.sqrt(3.0 / 5.0)
        t_c2j[7, 7] = np.sqrt(2.0 / 5.0)
        t_c2j[6, 8] = np.sqrt(4.0 / 5.0)
        t_c2j[9, 8] = np.sqrt(1.0 / 5.0)
        t_c2j[8, 9] = 1.0

        return t_c2j

    elif orb_l == 3:
        t_c2j = np.zeros((14, 14), dtype=np.complex128)
        t_c2j[0, 0] = -np.sqrt(6.0 / 7.0)
        t_c2j[3, 0] = np.sqrt(1.0 / 7.0)
        t_c2j[2, 1] = -np.sqrt(5.0 / 7.0)
        t_c2j[5, 1] = np.sqrt(2.0 / 7.0)
        t_c2j[4, 2] = -np.sqrt(4.0 / 7.0)
        t_c2j[7, 2] = np.sqrt(3.0 / 7.0)
        t_c2j[6, 3] = -np.sqrt(3.0 / 7.0)
        t_c2j[9, 3] = np.sqrt(4.0 / 7.0)
        t_c2j[8, 4] = -np.sqrt(2.0 / 7.0)
        t_c2j[11, 4] = np.sqrt(5.0 / 7.0)
        t_c2j[10, 5] = -np.sqrt(1.0 / 7.0)
        t_c2j[13, 5] = np.sqrt(6.0 / 7.0)
        t_c2j[1, 6] = 1.0
        t_c2j[0, 7] = np.sqrt(1.0 / 7.0)
        t_c2j[3, 7] = np.sqrt(6.0 / 7.0)
        t_c2j[2, 8] = np.sqrt(2.0 / 7.0)
        t_c2j[5, 8] = np.sqrt(5.0 / 7.0)
        t_c2j[4, 9] = np.sqrt(3.0 / 7.0)
        t_c2j[7, 9] = np.sqrt(4.0 / 7.0)
        t_c2j[6, 10] = np.sqrt(4.0 / 7.0)
        t_c2j[9, 10] = np.sqrt(3.0 / 7.0)
        t_c2j[8, 11] = np.sqrt(5.0 / 7.0)
        t_c2j[11, 11] = np.sqrt(2.0 / 7.0)
        t_c2j[10, 12] = np.sqrt(6.0 / 7.0)
        t_c2j[13, 12] = np.sqrt(1.0 / 7.0)
        t_c2j[12, 13] = 1.0

        return t_c2j

    else:
        raise Exception("error in tmat_c2j: Have NOT implemented for this case: ", orb_l)


def tmat_r2trig(case, ispin=False):
    """
    Get the unitary transformation matrix from the real spherical harmonics to the
    trigonal harmonic basis. See PHYSICAL REVIEW B 105, 245153 (2022) for the
    definition of the orbitals.

    Parameters
    ----------
    case: string
        Label for different systems.

        - '111': :math:`z` axis along diagonal of coordinate system
        - '100': :math:`z` axis vertical

    ispin: logical
        Whether to include spin degree of freedom or not (default: False).

    Returns
    -------
    t_r2trig: 2d complex array
        The transformation matrix.
    """
    t_r2trig = np.zeros((5, 5), dtype=complex)
    cone = np.complex128(1.0 + 0.0j)

    # trig-orbitals a1g epi_g epi_g esigma_g esigma_g
    # d-orbitals dz2, dzx, dzy, dx2-y2, dxy

    if case == '111':

        # a1g = 3z^2-r^2
        t_r2trig[0, 0] = cone

        # epi_g =  2 dxy/ sqrt(6) + dyz/sqrt(3)
        t_r2trig[4, 1] = 2*cone/np.sqrt(6)
        t_r2trig[2, 1] = cone/np.sqrt(3)

        # epi_g =  2 dx^2-y^2/ sqrt(6) - dzx/sqrt(3)
        t_r2trig[3, 2] = 2*cone/np.sqrt(6)
        t_r2trig[1, 2] = -cone/np.sqrt(3)

        # esigma_g =  dx^2-y^2/ sqrt(3) + 2dzx/sqrt(6)
        t_r2trig[3, 3] = cone/np.sqrt(3)
        t_r2trig[1, 3] = 2*cone/np.sqrt(6)

        # esigma_g =  dxy/ sqrt(3) - 2dyz/sqrt(6)
        t_r2trig[4, 4] = cone/np.sqrt(3)
        t_r2trig[2, 4] = -2*cone/np.sqrt(6)

    elif case == '100':
        t_r2trig[1, 0] = cone/np.sqrt(3)
        t_r2trig[2, 0] = cone/np.sqrt(3)
        t_r2trig[4, 0] = cone/np.sqrt(3)

        om = np.exp(1j*2*np.pi/3)
        t_r2trig[1, 1] = -cone/np.sqrt(3)
        t_r2trig[2, 1] = -om/np.sqrt(3)
        t_r2trig[4, 1] = -om**2/np.sqrt(3)

        t_r2trig[1, 2] = cone/np.sqrt(3)
        t_r2trig[2, 2] = om**-1/np.sqrt(3)
        t_r2trig[4, 2] = om**-2/np.sqrt(3)

        t_r2trig[1, 3] = -cone/np.sqrt(2)
        t_r2trig[3, 3] = -1j/np.sqrt(2)

        t_r2trig[1, 4] = cone/np.sqrt(2)
        t_r2trig[3, 4] = -1j/np.sqrt(2)
        
        # # a1g = (dzx + dzy + dxy)/sqrt(3)
        # t_r2trig[1, 0] = cone/np.sqrt(3)
        # t_r2trig[2, 0] = cone/np.sqrt(3)
        # t_r2trig[4, 0] = cone/np.sqrt(3)

        # # epi_g = (2dxy - dyz - dzx)/sqrt(6)
        # t_r2trig[4, 1] = 2/np.sqrt(6)
        # t_r2trig[2, 1] = -cone/np.sqrt(6)
        # t_r2trig[1, 1] = -cone/np.sqrt(6)

        # # epi_g = (dyz - dzx)/sqrt(2)
        # t_r2trig[2, 2] = cone/np.sqrt(2)
        # t_r2trig[1, 2] = -cone/np.sqrt(2)

        # # esigma_g =  dz^2
        # t_r2trig[0, 3] = cone

        # # esigma_g =  dx^2-y^2
        # t_r2trig[3, 4] = cone

    else:
        raise Exception(f"case {case} not implemented."
                        "It should be either '111' or '110'")

    # the spin order is: up dn up dn ... up dn
    if ispin:
        t_r2trig_spin = np.zeros((10, 10), dtype=complex)
        # spin up
        t_r2trig_spin[0::2, 0::2] = t_r2trig
        # spin dn
        t_r2trig_spin[1::2, 1::2] = t_r2trig
        return t_r2trig_spin
    else:
        return t_r2trig


def tmat_trig2r(case, ispin=False):
    """
    Get the unitary transformation matrix from the trigonal harmonic basis to the
    real spherical harmonics. See PHYSICAL REVIEW B 105, 245153 (2022) for the
    definition of the orbitals.

    Parameters
    ----------
    case: string
        Label for different systems.

        - '111': :math:`z` axis along diagonal of coordinate system
        - '100': :math:`z` axis vertical

    ispin: logical
        Whether to include spin degree of freedom or not (default: False).

    Returns
    -------
    t_trig2r: 2d complex array
        The transformation matrix.
    """
    t_trig2r = np.conj(np.transpose(tmat_r2trig(case, ispin=ispin)))
    return t_trig2r


def transform_utensor(umat, tmat):
    """
    Transform the rank-4 Coulomb interaction tensor from one basis to
    another basis.

    Parameters
    ----------
    umat: 4d array
        Coulomb interaction tensor (rank-4).
    tmat: 2d array
        The transformation matrix.

    Returns
    -------
    umat_new: 4d complex array
        The Coulomb interaction tensor in the new basis.
    """

    n = umat.shape[0]
    umat_new = np.zeros((n, n, n, n), dtype=np.complex128)

    a1, a2, a3, a4 = np.nonzero(abs(umat) > 1E-16)
    nonzero = np.stack((a1, a2, a3, a4), axis=-1)

    for ii, jj, kk, mm in nonzero:
        for i in range(n):
            if abs(tmat[ii, i]) < 1E-16:
                continue
            else:
                for j in range(n):
                    if abs(tmat[jj, j]) < 1E-16:
                        continue
                    else:
                        for k in range(n):
                            if abs(tmat[kk, k]) < 1E-16:
                                continue
                            else:
                                for m in range(n):
                                    umat_new[i, j, k, m] += (np.conj(tmat[ii, i]) *
                                                             np.conj(tmat[jj, j]) *
                                                             umat[ii, jj, kk, mm] *
                                                             tmat[kk, k] *
                                                             tmat[mm, m])

    return umat_new


def fourier_hr2hk(norbs, nkpt, kvec, nrpt, rvec, deg_rpt, hr):
    """
    Fourier transform a tight-binding Hamiltonian :math:`H(r)` from
    real space to :math:`k` space :math:`H(k)`,
    for Wannier90.

    Parameters
    ----------
    norbs: int
        Number of orbitals.
    nkpt: int
        Number of :math:`k`-points.
    kvec: 2d float array
        Fractional coordinate for k-points.
    nrpt: int
        Number of r-points.
    rvec: 2d float array
        Fractional coordinate for r-points.
    deg_rpt: int
        The degenerancy for each r-point.
    hr: 3d complex array
        Hamiltonian in r-space.

    Returns
    -------
    hk: 3d complex array
        Hamiltonian in k-space.
    """

    hk = np.zeros((nkpt, norbs, norbs), dtype=np.complex128)
    for i in range(nkpt):
        for j in range(nrpt):
            coef = -2 * np.pi * np.dot(kvec[i, :], rvec[j, :])
            ratio = (np.cos(coef) + np.sin(coef) * 1j) / float(deg_rpt[j])
            hk[i, :, :] = hk[i, :, :] + ratio * hr[j, :, :]
    return hk
