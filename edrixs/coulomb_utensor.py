__all__ = ['get_gaunt', 'umat_slater', 'get_umat_kanamori_ge', 'get_F0',
           'get_umat_slater', 'get_umat_kanamori', 'get_umat_slater_3shells']

import numpy as np
from sympy.physics.wigner import gaunt
from .basis_transform import tmat_c2r, tmat_r2c, tmat_c2j, transform_utensor
from .utils import info_atomic_shell, case_to_shell_name, slater_integrals_name


def get_gaunt(l1, l2):
    """
    Calculate the Gaunt coefficents :math:`C_{l_1,l_2}(k,m_1,m_2)`

    .. math::

        C_{l_1,l_2}(k,m_1,m_2)=\\sqrt{\\frac{4\\pi}{2k+1}} \\int
        \\mathop{d\\phi} \\mathop{d\\theta} sin(\\theta)
        Y_{l_1}^{m_1\\star}(\\theta,\\phi) Y_{k}^{m_1-m_2}(\\theta,\\phi)
        Y_{l_2}^{m_2}(\\theta,\\phi)

    Parameters
    ----------
    l1: int
        The first quantum number of angular momentum.
    l2: int
        The second quantum number of angular momentum.

    Returns
    -------
    res: 3d float array
        The calculated Gaunt coefficents.

        The 1st index (:math:`= 0, 1, ..., l_1+l_2+1`) is the order :math:`k`.

        The 2nd index (:math:`= 0, 1, ... ,2l_1`) is the magnetic quantum
        number :math:`m_1` plus :math:`l_1`

        The 3nd index (:math:`= 0, 1, ... ,2l_2`) is the magnetic quantum
        number :math:`m_2` plus :math:`l_2`

    Notes
    -----
    It should be noted that :math:`C_{l_1,l_2}(k,m_1,m_2)` is
    nonvanishing only when

    :math:`k + l_1 + l_2 = \\text{even}`,

    and

    :math:`|l_1 -  l_2| \\leq k \\leq l_1 + l_2`.

    Please see Ref. [1]_ p. 10 for more details.

    References
    ----------
    .. [1] Sugano S, Tanabe Y and Kamimura H. 1970. Multiplets of
       Transition-Metal Ions in Crystals. Academic Press, New York and London.

    Examples
    --------
    >>> import edrixs

    Get gaunt coefficients between :math:`p`-shell and :math:`d`-shell

    >>> g = edrixs.get_gaunt(1, 2)

    """

    from sympy import N
    res = np.zeros((l1 + l2 + 1, 2 * l1 + 1, 2 * l2 + 1), dtype=np.float64)
    for k in range(l1 + l2 + 1):
        if not (np.mod(l1 + l2 + k, 2) == 0 and np.abs(l1 - l2) <= k <= l1 + l2):
            continue
        for i1, m1 in enumerate(range(-l1, l1 + 1)):
            for i2, m2 in enumerate(range(-l2, l2 + 1)):
                res[k, i1, i2] = (N(gaunt(l1, k, l2, -m1, m1 - m2, m2)) *
                                  (-1.0)**m1 * np.sqrt(4 * np.pi / (2 * k + 1)))
    return res


def umat_slater(l_list, fk):
    """
    Calculate the Coulomb interaction tensor
    
    .. math::

        \\lange i j | U_{m_{l_i}m_{s_i}, m_{l_j}m_{s_j}, m_{l_t}m_{s_t},
        m_{l_u}m_{s_u}}^{i,j,t,u} |t u \\rangle

    which is parameterized by
    Slater integrals :math:`F^{k}`:

    .. math::

        U_{m_{l_i}m_{s_i}, m_{l_j}m_{s_j}, m_{l_t}m_{s_t},
        m_{l_u}m_{s_u}}^{i,j,t,u}
        =\\frac{1}{2} \\delta_{m_{s_i},m_{s_t}}\\delta_{m_{s_j},m_{s_u}}
        \\delta_{m_{l_i}+m_{l_j}, m_{l_t}+m_{l_u}}
        \\sum_{k}C_{l_i,l_t}(k,m_{l_i},m_{l_t})C_{l_u,l_j}
        (k,m_{l_u},m_{l_j})F^{k}_{i,j,t,u}

    where :math:`m_s` is the magnetic quantum number for spin
    and :math:`m_l` is the magnetic quantum number for orbital.
    :math:`F^{k}_{i,j,t,u}` are Slater integrals.
    :math:`C_{l_i,l_j}(k,m_{l_i},m_{l_j})` are Gaunt coefficients.

    Parameters
    ----------
    l_list: list of int
        contains the quantum number of orbital angular momentum
        :math:`l` for each shell.
    fk: dict of float
        contains all the possible Slater integrals between the shells
        in l_list, the key is a tuple of 5 ints (:math:`k,i,j,t,u`),
        where :math:`k` is the order,
        :math:`i,j,t,u` are the shell indices begin with 1.

    Returns
    -------
    umat: 4d array of complex
        contains the Coulomb interaction tensor.

    Examples
    --------
    >>> import edrixs

    For only one :math:`d`-shell

    >>> l_list = [2]
    >>> fk={}
    >>> F0, F2, F4 = 5.0, 4.0 2.0
    >>> fk[(0,1,1,1,1)] = F0
    >>> fk[(2,1,1,1,1)] = F2
    >>> fk[(4,1,1,1,1)] = F4
    >>> umat_d = edrixs.umat_slater(l_list, fk)

    For one :math:`d`-shell and one :math:`p`-shell

    >>> l_list = [2,1]
    >>> fk={}
    >>> F0_dd, F2_dd, F4_dd = 5.0, 4.0, 2.0
    >>> F0_dp, F2_dp = 4.0, 2.0
    >>> G1_dp, G3_dp = 2.0, 1.0
    >>> F0_pp, F2_pp = 2.0, 1.0

    >>> fk[(0,1,1,1,1)] = F0_dd
    >>> fk[(2,1,1,1,1)] = F2_dd
    >>> fk[(4,1,1,1,1)] = F4_dd

    >>> fk[(0,1,2,1,2)] = F0_dp
    >>> fk[(0,2,1,2,1)] = F0_dp
    >>> fk[(2,1,2,1,2)] = F2_dp
    >>> fk[(2,2,1,2,1)] = F2_dp

    >>> fk[(1,1,2,2,1)] = G1_dp
    >>> fk[(1,2,1,1,2)] = G1_dp
    >>> fk[(3,1,2,2,1)] = G3_dp
    >>> fk[(3,2,1,1,2)] = G3_dp

    >>> fk[(0,2,2,2,2)] = F0_pp
    >>> fk[(2,2,2,2,2)] = F2_pp

    >>> umat_dp = edrixs.umat_slater(l_list, fk)

    See also
    --------
    coulomb_utensor.get_umat_slater
    coulomb_utensor.get_umat_kanamori
    coulomb_utensor.get_umat_kanamori_ge
    """
    k_list = list(range(0, 2 * max(l_list) + 1))
    ck = {}
    orb_label = []
    # for each shell
    for i, l in enumerate(l_list):
        # magnetic quantum number
        for m in range(-l, l + 1):
            # spin part, up,dn
            for spin in range(2):
                orb_label.append((i + 1, l, m, spin))

    # all the possible gaunt coefficient
    for l1 in l_list:
        for l2 in l_list:
            ck_tmp = get_gaunt(l1, l2)
            for m1 in range(-l1, l1 + 1):
                for m2 in range(-l2, l2 + 1):
                    for k in k_list:
                        if (np.mod(l1 + l2 + k, 2) == 0 and np.abs(l1 - l2) <= k <= l1 + l2):
                            ck[(k, l1, m1, l2, m2)] = ck_tmp[k, m1 + l1, m2 + l2]
                        else:
                            ck[(k, l1, m1, l2, m2)] = 0

    # build the coulomb interaction tensor
    norbs = len(orb_label)
    umat = np.zeros((norbs, norbs, norbs, norbs), dtype=np.complex128)
    for orb1 in range(norbs):
        for orb2 in range(norbs):
            for orb3 in range(norbs):
                for orb4 in range(norbs):
                    i1, l1, m1, sigma1 = orb_label[orb1]
                    i2, l2, m2, sigma2 = orb_label[orb2]
                    i3, l3, m3, sigma3 = orb_label[orb3]
                    i4, l4, m4, sigma4 = orb_label[orb4]
                    if (m1 + m2) != (m3 + m4):
                        continue
                    if (sigma1 != sigma3) or (sigma2 != sigma4):
                        continue
                    res = 0.0
                    for k in k_list:
                        tmp_key = (k, i1, i2, i3, i4)
                        if tmp_key in list(fk.keys()):
                            res += ck[(k, l1, m1, l3, m3)] * ck[(k, l4, m4, l2, m2)] * fk[tmp_key]
                    umat[orb1, orb2, orb4, orb3] = res

    umat = umat / 2.0
    return umat


def get_umat_kanamori_ge(norbs, U1, U2, J, Jx, Jp):
    """
    Calculate the Coulomb interaction tensor for a Kanamori-type interaction.
    For the general case, it is parameterized by :math:`U_1, U_2, J, J_x, J_p`.

    Parameters
    ----------
    norbs: int
        number of orbitals (including spin).
    U1: float
        Hubbard :math:`U` for electrons residing on the same orbital with opposite spin.
    U2: float
        Hubbard :math:`U` for electrons residing on different orbitals.
    J: float
        Hund's coupling for density-density interaction.
    Jx: float
        Hund's coupling for spin flip.
    Jp: float
        Hund's coupling for pair-hopping.

    Returns
    -------
    umat: 4d complex array
        The calculated Coulomb interaction tensor

    Notes
    -----
    The order of spin index is: up, down, up, down, ..., up, down.

    See Also
    --------
    coulomb_utensor.get_umat_kanamori
    coulomb_utensor.get_umat_slater
    coulomb_utensor.umat_slater

    """

    umat = np.zeros((norbs, norbs, norbs, norbs), dtype=np.complex128)
    for alpha in range(0, norbs - 1):
        for beta in range(alpha + 1, norbs):
            for gamma in range(0, norbs - 1):
                for delta in range(gamma + 1, norbs):
                    aband = alpha // 2
                    aspin = alpha % 2
                    bband = beta // 2
                    bspin = beta % 2
                    gband = gamma // 2
                    gspin = gamma % 2
                    dband = delta // 2
                    dspin = delta % 2

                    dtmp = 0.0
                    if ((alpha == gamma) and (beta == delta)):
                        if ((aband == bband) and (aspin != bspin)):
                            dtmp = dtmp + U1

                    if ((alpha == gamma) and (beta == delta)):
                        if (aband != bband):
                            dtmp = dtmp + U2

                    if ((alpha == gamma) and (beta == delta)):
                        if ((aband != bband) and (aspin == bspin)):
                            dtmp = dtmp - J

                    if ((aband == gband) and (bband == dband)):
                        if ((aspin != gspin) and (bspin != dspin)
                                and (aspin != bspin)):
                            dtmp = dtmp - Jx

                    if ((aband == bband) and (dband == gband)
                            and (aband != dband)):
                        if ((aspin != bspin) and (dspin != gspin)
                                and (aspin == gspin)):
                            dtmp = dtmp + Jp

                    umat[alpha, beta, delta, gamma] = dtmp

    return umat


def get_F0(case, *args):
    case = case.strip()
    if case == 's':
        return 0.0

    elif case == 'p':
        F2 = args[0]
        return 2.0 / 25.0 * F2

    elif case == 'd':
        F2, F4 = args[0:2]
        return 2.0 / 63.0 * F2 + 2.0 / 63.0 * F4

    elif case == 'f':
        F2, F4, F6 = args[0:3]
        return 4.0 / 195.0 * F2 + 2.0 / 143.0 * F4 + 100.0 / 5577.0 * F6

    elif case == 'ss':
        G0 = args[0]
        return 0.5 * G0

    elif case == 'sp' or case == 'ps':
        G1 = args[0]
        return 1.0 / 6.0 * G1

    elif case == 'sd' or case == 'ds':
        G2 = args[0]
        return 1.0 / 10.0 * G2

    elif case == 'sf' or case == 'fs':
        G3 = args[0]
        return 1.0 / 14.0 * G3

    elif case == 'pp':
        G0, G2 = args[0:2]
        return 1.0 / 6.0 * G0 + 1.0 / 15.0 * G2

    elif case == 'pd' or case == 'dp':
        G1, G3 = args[0:2]
        return 1.0 / 15.0 * G1 + 3.0 / 70.0 * G3

    elif case == 'pf' or case == 'fp':
        G2, G4 = args[0:2]
        return 3.0 / 70.0 * G2 + 2.0 / 63.0 * G4

    elif case == 'dd':
        G0, G2, G4 = args[0:3]
        return 1.0 / 10.0 * G0 + 1.0 / 35.0 * G2 + 1.0 / 35.0 * G4

    elif case == 'df' or case == 'fd':
        G1, G3, G5 = args[0:3]
        return 3.0 / 70.0 * G1 + 2.0 / 105.0 * G3 + 5.0 / 231.0 * G5

    elif case == 'ff':
        G0, G2, G4, G6 = args[0:4]
        return 1.0 / 14.0 * G0 + 2.0 / 105.0 * G2 + 1.0 / 77.0 * G4 + 50.0 / 3003.0 * G6

    else:
        raise Exception("error in get_F0():  Unknown case name:", case)


def get_umat_slater(case, *args):
    """
    Convenient adapter function to return the Coulomb interaction tensor for common case.

    Parameters
    ----------
    case: string
        Indicates atomic shells, should be one of

        For single shell:

        - 's':    single :math:`s`-shell (:math:`l=0`)
        - 'p':    single :math:`p`-shell (:math:`l=1`)
        - 'p12':  single :math:`p_{1/2}`-shell (:math:`l=1`)
        - 'p32':  single :math:`p_{3/2}`-shell (:math:`l=1`)
        - 't2g':  single :math:`t_{2g}`-shell (:math:`l_{\\text{eff}}=1`)
        - 'd':    single :math:`d`-shell (:math:`l=2`)
        - 'd32':  single :math:`d_{3/2}`-shell (:math:`l=2`)
        - 'd52':  single :math:`d_{5/2}`-shell (:math:`l=2`)
        - 'f':    single :math:`f`-shell (:math:`l=3`)
        - 'f52':  single :math:`f_{5/2}`-shell (:math:`l=3`)
        - 'f72':  single :math:`f_{7/2}`-shell (:math:`l=3`)

        For two shells:

        case = str1 + str2

        where, str1 and str2 are strings and they can be any of

        ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']

        For examples,

        - 'dp':     1st :math:`d`-shell and 2nd :math:`p`-shell
        - 'dp12':   1st :math:`d`-shell and 2nd :math:`p_{1/2}`-shell
        - 'f52p32': 1st :math:`f_{5/2}`-shell and 2nd :math:`p_{3/2}`-shell
        - 't2gp':   1st :math:`t_{2g}`-shell and 2nd :math:`p`-shell

    *args: floats
        Variable length argument list. Slater integrals.
        The order of these integrals shoule be

        For only one shell case,

        args = [F0, F2, F4, F6, ....]

        For two shells case,

        args = [FX_11, FX_12, GX_12, FX_22]

        where, 1 (2) means 1st (2nd)-shell, and X=0, 2, 4, ... or X=1, 3, 5 ...,
        and X should be in ascending order. The following are possible cases:

        - 's':

        args = [F0]

        - 'p', 'p12', 'p32':

        args = [F0, F2]

        - 'd', 'd32', 'd52', 't2g':

        args = [F0, F2, F4]

        - 'f', 'f52', 'f72':

        args = [F0, F2, F4, F6]

        - 'ss':

        args = [F0_11, F0_12, G0_12, F0_22]

        - 'ps', 'p12s', 'p32s':

        args = [F0_11, F2_11, F0_12, G1_12, F0_22]

        - 'ds', 'd32s', 'd52s', 't2gs':

        args = [F0_11, F2_11, F4_11, F0_12, G2_12, F0_22]

        - 'fs', 'f52s', 'f72s':

        args = [F0_11, F2_11, F4_11, F6_11, F0_12, G3_12, F0_22]

        - 'sp', 'sp12', 'sp32':

        args = [F0_11, F0_12, G1_12, F0_22, F2_22]

        - 'pp', 'pp12', 'pp32', 'p12p', 'p12p12', 'p12p32', 'p32p', 'p32p12', 'p32p32':

        args = [F0_11, F2_11, F0_12, F2_12, G0_12, G2_12, F0_22, F2_22]

        - 'dp', 'dp12', 'dp32', 'd32p', 'd32p12', 'd32p32', 'd52p', 'd52p12', 'd52p32',
          't2gp', 't2gp12', t2gp32':

        args = [F0_11, F2_11, F4_11, F0_12, F2_12, G1_12, G3_12, F0_22, F2_22]

        - 'fp', 'fp12', 'fp32', 'f52p', 'f52p12', 'f52p32', 'f72p', 'f72p12', 'f72p32':

        args = [F0_11, F2_11, F4_11, F6_11, F0_12, F2_12, G2_12, G4_12, F0_22, F2_22]

        - 'sd', 'sd32', 'sd52':

        args = [F0_11, F0_12, G2_12, F0_22, F2_22, F4_22]

        - 'pd', 'pd32', 'pd52', 'p12d', 'p12d32', 'p12d52', 'p32d', 'p32d32', 'p32d52':

        args = [F0_11, F2_11, F0_12, F2_12, G1_12, G3_12, F0_22, F2_22, F4_22]

        - 'dd', 'dd32', 'dd52', 'd32d', 'd32d32', 'd32d52', 'd52d', 'd52d32', 'd52d52',
          't2gd', 't2gd32', 't2gd52':

        args = [F0_11, F2_11, F4_11, F0_12, F2_12, F4_12, G0_12, G2_12, G4_12, F0_22, F2_22, F4_22]

        - 'fd', 'fd32', 'fd52', 'f52d', 'f52d32', 'f52d52', 'f72d', 'f72d32', 'f72d52':

        args = [F0_11, F2_11, F4_11, F6_11, F0_12, F2_12, F4_12, G1_12, G3_12, G5_12,
        F0_22, F2_22, F4_22]

        - 'sf', 'sf52', 'sf72':

        args = [F0_11, F0_12, G3_12, F0_22, F2_22, F4_22, F6_22]

        - 'pf', 'pf52', 'pf72', 'p12f', 'p12f52', 'p12f72', 'p32f', 'p32f52', 'p32f72':

        args = [F0_11, F2_11, F0_12, F2_12, G2_12, G4_12, F0_22, F2_22, F4_22, F6_22]

        - 'df', 'df52', 'df72', 'd32f', 'd32f52', 'd32f72', 'd52f', 'd52f52', 'd52f72',
          't2gf', 't2gf52', 't2gf72':

        args = [F0_11, F2_11, F4_11, F0_12, F2_12, F4_12, G1_12, G3_12, G5_12,
        F0_22, F2_22, F4_22, F6_22]

        - 'ff', 'ff52', 'ff72', 'f52f', 'f52f52', 'f52f72', 'f72f', 'f72f52', 'f72f72':

        args = [F0_11, F2_11, F4_11, F6_11, F0_12, F2_12, F4_12, F6_12,
        G0_12, G2_12, G4_12, G6_12, F0_22, F2_22, F4_22, F6_22]

    Returns
    -------
    umat: 4d array of complex
        the Coulomb interaction tensor

    Examples
    --------
    >>> import edrixs
    >>> F0_dd, F2_dd, F4_dd = 3.0, 1.0, 0.5
    >>> F0_dp, F2_dp, G1_dp, G3_dp = 2.0, 1.0, 0.2, 0.1
    >>> F0_pp, F2_pp = 0.0, 0.0
    >>> slater = [F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp, F0_pp, F2_pp]
    >>> umat_d = edrixs.get_umat_slater('d', F0_dd, F2_dd, F4_dd)
    >>> umat_dp = edrixs.get_umat_slater('dp', *slater)
    >>> umat_t2gp = edrixs.get_umat_slater('t2gp', *slater)
    >>> umat_dp32 = edrixs.get_umat_slater('dp32', *slater)

    See Also
    --------
    coulomb_utensor.umat_slater
    coulomb_utensor.get_umat_kanamori
    coulomb_utensor.get_umat_kanamori_ge

    """
    info = info_atomic_shell()
    shells = case_to_shell_name(case)
    nslat = len(slater_integrals_name(shells))
    if nslat != len(args):
        raise Exception("Number of Slater integrals", len(args), " is not equal to ", nslat)

    special_shell = ['t2g', 'p12', 'p32', 'd32', 'd52', 'f52', 'f72']
    orb_indx = {
        special_shell[0]: [2, 3, 4, 5, 8, 9],
        special_shell[1]: [0, 1],
        special_shell[2]: [2, 3, 4, 5],
        special_shell[3]: [0, 1, 2, 3],
        special_shell[4]: [4, 5, 6, 7, 8, 9],
        special_shell[5]: [0, 1, 2, 3, 4, 5],
        special_shell[6]: [6, 7, 8, 9, 10, 11, 12, 13]
    }
    # only one shell
    if len(shells) == 1:
        orbl = info[shells[0]][0]
        l_list = [orbl]
        fk = {}
        it = 0
        for rank in range(0, 2*orbl+1, 2):
            fk[(rank, 1, 1, 1, 1)] = args[it]
            it += 1
        umat = umat_slater(l_list, fk)
        # truncate to a sub-shell if necessary
        if shells[0] in special_shell:
            if shells[0] == 't2g':
                tmat = tmat_c2r('d', True)
            else:
                tmat = tmat_c2j(orbl)
            umat = transform_utensor(umat, tmat)
            indx = orb_indx[shells[0]]
            umat_tmp = np.zeros((len(indx), len(indx), len(indx), len(indx)), dtype=complex)
            umat_tmp[:, :, :, :] = umat[indx][:, indx][:, :, indx][:, :, :, indx]
            if shells[0] == 't2g':
                umat_tmp[:, :, :, :] = transform_utensor(umat_tmp, tmat_r2c('t2g', True))
            umat = umat_tmp
    # two shells
    elif len(shells) == 2:
        name1, name2 = shells
        l1, l2 = info[name1][0], info[name2][0]
        n1, n2 = 2*(2*l1+1), 2*(2*l2+1)
        ntot = n1 + n2
        l_list = [l1, l2]
        fk = {}
        it = 0
        for rank in range(0, 2 * l1 + 1, 2):
            fk[(rank, 1, 1, 1, 1)] = args[it]
            it += 1
        for rank in range(0, min(2 * l1, 2 * l2) + 1, 2):
            fk[(rank, 1, 2, 1, 2)] = args[it]
            fk[(rank, 2, 1, 2, 1)] = args[it]
            it += 1
        for rank in range(abs(l1 - l2), l1 + l2 + 1, 2):
            fk[(rank, 1, 2, 2, 1)] = args[it]
            fk[(rank, 2, 1, 1, 2)] = args[it]
            it += 1
        for rank in range(0, 2 * l2 + 1, 2):
            fk[(rank, 2, 2, 2, 2)] = args[it]
            it += 1
        umat = umat_slater(l_list, fk)
        # truncate to a sub-shell if necessary
        if (name1 in special_shell) or (name2 in special_shell):
            tmat = np.eye(ntot, dtype=complex)
            indx1 = list(range(0, n1))
            if name1 in special_shell:
                if name1 == 't2g':
                    tmat[0:n1, 0:n1] = tmat_c2r('d', True)
                else:
                    tmat[0:n1, 0:n1] = tmat_c2j(l1)
                indx1 = orb_indx[name1]

            indx2 = [n1 + i for i in range(0, n2)]
            if name2 in special_shell:
                if name2 == 't2g':
                    tmat[n1:ntot, n1:ntot] = tmat_c2r('d', True)
                else:
                    tmat[n1:ntot, n1:ntot] = tmat_c2j(l2)
                indx2 = [n1 + i for i in orb_indx[name2]]

            indx = indx1 + indx2
            umat = transform_utensor(umat, tmat)
            umat_tmp = np.zeros((len(indx), len(indx), len(indx), len(indx)), dtype=complex)
            umat_tmp[:, :, :, :] = umat[indx][:, indx][:, :, indx][:, :, :, indx]
            if name1 == 't2g' or name2 == 't2g':
                tmat = np.eye(len(indx), dtype=np.complex128)
                if name1 == 't2g':
                    tmat[0:6, 0:6] = tmat_r2c('t2g', True)
                if name2 == 't2g':
                    tmat[-6:, -6:] = tmat_r2c('t2g', True)
                umat_tmp[:, :, :, :] = transform_utensor(umat_tmp, tmat)
            umat = umat_tmp
    else:
        raise Exception("Not implemented for this case: ", shells)

    return umat


def get_umat_slater_3shells(shell_name, *args):
    """
    Given three shells, build the slater type of Coulomb tensors among
    the three shells.

    Parameters
    ----------
    shell_name: tuple of three strings
        Shells names.
    *args: floats
        Slater integrals. The order should be

        FX_11, FX_12, GX_12, FX_22, FX_13, GX_13, FX_23, GX_23, FX_33

        where, 1, 2, 3 means 1st, 2nd, 3rd shell, and X=0, 2, 4, ... or X=1, 3, 5 ...,
        and X should be in ascending order.

    Returns
    -------
    umat: 4d complex array
        Rank-4 Coulomb tensors.
    """

    v1_name = shell_name[0].strip()
    v2_name = shell_name[1].strip()
    v3_name = shell_name[2].strip()
    info_shell = info_atomic_shell()

    v1_orbl = info_shell[v1_name][0]
    v2_orbl = info_shell[v2_name][0]
    v3_orbl = info_shell[v3_name][0]

    v1_norb = info_shell[v1_name][1]
    v2_norb = info_shell[v2_name][1]
    v3_norb = info_shell[v3_name][1]

    # total number of orbitals
    ntot = v1_norb + v2_norb + v3_norb
    v1v2_norb = v1_norb + v2_norb

    indx = []
    it = 0
    it += len(range(0, 2*v1_orbl+1, 2))
    indx.append(it)
    it += len(range(0, min(2*v1_orbl, 2*v2_orbl)+1, 2))
    it += len(range(abs(v1_orbl-v2_orbl), v1_orbl+v2_orbl+1, 2))
    indx.append(it)
    it += len(range(0, 2*v2_orbl+1, 2))
    indx.append(it)
    it += len(range(0, min(2*v1_orbl, 2*v3_orbl)+1, 2))
    it += len(range(abs(v1_orbl-v3_orbl), v1_orbl+v3_orbl+1, 2))
    indx.append(it)
    it += len(range(0, min(2*v2_orbl, 2*v3_orbl)+1, 2))
    it += len(range(abs(v2_orbl-v3_orbl), v2_orbl+v3_orbl+1, 2))
    indx.append(it)
    it += len(range(0, 2*v3_orbl+1, 2))
    indx.append(it)

    if it != len(args):
        raise Exception("Number of Slater integrals", len(args), " is not equal to ", it)

    umat = np.zeros((ntot, ntot, ntot, ntot), dtype=complex)

    # v1-v2
    case = v1_name + v2_name
    arg_list = list(args[0:indx[0]]) + list(args[indx[0]:indx[1]]) + list(args[indx[1]:indx[2]])
    umat_tmp = get_umat_slater(case, *arg_list)
    umat[0:v1v2_norb, 0:v1v2_norb, 0:v1v2_norb, 0:v1v2_norb] = umat_tmp

    # v1-v3
    case = v1_name + v3_name
    arg_list = [0.0] * indx[0] + list(args[indx[2]:indx[3]]) + list(args[indx[4]:indx[5]])
    umat_tmp = get_umat_slater(case, *arg_list)
    aa = list(range(0, v1_norb)) + list(range(v1v2_norb, ntot))
    for i in range(v1_norb + v3_norb):
        for j in range(v1_norb + v3_norb):
            for k in range(v1_norb + v3_norb):
                for m in range(v1_norb + v3_norb):
                    umat[aa[i], aa[j], aa[k], aa[m]] += umat_tmp[i, j, k, m]

    # v2-v3
    case = v2_name + v3_name
    arg_list = ([0.0] * (indx[2] - indx[1]) +
                list(args[indx[3]:indx[4]]) +
                [0.0] * (indx[5] - indx[4]))
    umat_tmp = get_umat_slater(case, *arg_list)
    aa = list(range(v1_norb, ntot))
    for i in range(v2_norb + v3_norb):
        for j in range(v2_norb + v3_norb):
            for k in range(v2_norb + v3_norb):
                for m in range(v2_norb + v3_norb):
                    umat[aa[i], aa[j], aa[k], aa[m]] += umat_tmp[i, j, k, m]

    return umat


def get_umat_kanamori(norbs, U, J):
    """
    Calculate the Coulomb interaction tensor for a Kanamori-type interaction.
    For the :math:`t2g`-shell case, it is parameterized by :math:`U, J`.

    Parameters
    ----------
    norbs: int
        number of orbitals (including spin).
    U: float
        Hubbard :math:`U` for electrons residing on the same orbital with opposite spin.
    J: float
        Hund's coupling.

    Returns
    -------
    umat: 4d complex array
        The calculated Coulomb interaction tensor

    Notes
    -----
    The order of spin index is: up, down, up, down, ..., up, down.

    See Also
    --------
    coulomb_utensor.get_umat_kanamori_ge
    coulomb_utensor.get_umat_slater
    coulomb_utensor.umat_slater

    """

    return get_umat_kanamori_ge(norbs, U, U - 2 * J, J, J, J)
