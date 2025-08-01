__all__ = ['beta_to_kelvin', 'kelvin_to_beta', 'boltz_dist', 'UJ_to_UdJH',
           'UdJH_to_UJ', 'UdJH_to_F0F2F4', 'UdJH_to_F0F2F4F6', 'F0F2F4_to_UdJH',
           'F0F2F4_to_UJ', 'F0F2F4F6_to_UdJH', 'CT_imp_bath',
           'CT_imp_bath_core_hole', 'info_atomic_shell',
           'case_to_shell_name', 'edge_to_shell_name', 'slater_integrals_name',
           'get_atom_data', 'rescale']

import numpy as np
import json
from importlib import resources


def beta_to_kelvin(beta):
    """
    Convert :math:`\\beta` to Kelvin.

    Parameters
    ----------
    beta: float
        Inversion temperature.

    Returns
    -------
    T: float
        Temperature (K).
    """

    kb = 8.6173303E-5
    ev = 1.0 / float(beta)
    T = ev / kb
    return T


def kelvin_to_beta(k):
    """
    Convert temperature from Kelvin to :math:`\\beta`.

    Parameters
    ----------
    k: float
        Temperature in Kelvin.

    Returns
    -------
    beta: float
        Inversion temperature.
    """

    kb = 8.6173303E-5
    beta = 1.0 / (kb * k)
    return beta


def boltz_dist(gs, T):
    """
    Return Boltzmann distributition.

    Parameters
    ----------
    gs: 1d float array
        Energy levels.
    T: float
        Temperature in Kelvin.

    Returns
    -------
    res: 1d float array
        The Boltzmann distributition.
    """
    tmp_gs = np.array(gs)
    beta = kelvin_to_beta(T)
    res = np.exp(-beta * (tmp_gs - min(tmp_gs))) / np.sum(np.exp(-beta * (tmp_gs - min(tmp_gs))))
    return res


def UJ_to_UdJH(U, J):
    """
    Given Kanamori :math:`U` and :math:`J`, return :math:`U_d` and :math:`J_H`,
    for :math:`t2g`-orbitals.

    Parameters
    ----------
    U: float
        Coulomb interaction :math:`U`.
    J: float
        Hund's coupling :math:`J`.

    Returns
    -------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_{H}`.
    """

    F2 = J / (3.0 / 49.0 + 20 * 0.625 / 441.0)
    F4 = 0.625 * F2
    JH = (F2 + F4) / 14.0
    Ud = U - 4.0 / 49.0 * (F2 + F4)

    return Ud, JH


def UdJH_to_UJ(Ud, JH):
    """
    Given :math:`U_d` and :math:`J_H`, return Kanamori :math:`U` and :math:`J`,
    for :math:`t2g`-orbitals.

    Parameters
    ----------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_H`.

    Returns
    -------
    U: float
        Coulomb interaction :math:`U` in Kanamori form.
    J: float
        Hund's coupling :math:`J` in Kanamori form.
    """

    F2 = 14.0 / 1.625 * JH
    F4 = 0.625 * F2
    J = 3.0 / 49.0 * F2 + 20 / 441.0 * F4
    U = Ud + 4.0 / 49.0 * (F2 + F4)

    return U, J


def UdJH_to_F0F2F4(Ud, JH):
    """
    Given :math:`U_d` and :math:`J_H`, return :math:`F_0`, :math:`F_2` and :math:`F_4`,
    for :math:`d`-orbitals.

    Parameters
    ----------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_H`.

    Returns
    -------
    F0: float
        Slater integral :math:`F_0`.
    F2: float
        Slater integral :math:`F_2`.
    F4: float
        Slater integral :math:`F_4`.
    """

    F0 = Ud
    F2 = 14 / 1.625 * JH
    F4 = 0.625 * F2

    return F0, F2, F4


def UdJH_to_F0F2F4F6(Ud, JH):
    """
    Given :math:`U_d` and :math:`J_H`, return :math:`F_0`, :math:`F_2`, :math:`F_4` and :math:`F_6`,
    for :math:`f`-orbitals.

    Parameters
    ----------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_H`.

    Returns
    -------
    F0: float
        Slater integral :math:`F_0`.
    F2: float
        Slater integral :math:`F_2`.
    F4: float
        Slater integral :math:`F_4`.
    F6: float
        Slater integral :math:`F_6`.

    """

    F0 = Ud
    F2 = 6435 / (286.0 + (195 * 451) / 675.0 + (250 * 1001) / 2025.0) * JH
    F4 = 451 / 675.0 * F2
    F6 = 1001 / 2025.0 * F2

    return F0, F2, F4, F6


def F0F2F4_to_UdJH(F0, F2, F4):
    """
    Given :math:`F_0`, :math:`F_2` and :math:`F_4`, return :math:`U_d` and :math:`J_H`,
    for :math:`d`-orbitals.

    Parameters
    ----------
    F0: float
        Slater integral :math:`F_0`.
    F2: float
        Slater integral :math:`F_2`.
    F4: float
        Slater integral :math:`F_4`.

    Returns
    -------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_H`.
    """

    Ud = F0
    JH = (F2 + F4) / 14.0

    return Ud, JH


def F0F2F4_to_UJ(F0, F2, F4):
    """
    Given :math:`F_0`, :math:`F_2` and :math:`F_4`, return :math:`U` and :math:`J`,
    for :math:`t2g`-orbitals.

    Parameters
    ----------
    F0: float
        Slater integral :math:`F_0`.
    F2: float
        Slater integral :math:`F_2`.
    F4: float
        Slater integral :math:`F_4`.

    Returns
    -------
    U: float
        Coulomb interaction :math:`U`.
    J: float
        Hund's coupling :math:`J`.
    """

    U = F0 + 4.0 / 49.0 * (F2 + F4)
    J = 3.0 / 49.0 * F2 + 20 / 441.0 * F4

    return U, J


def F0F2F4F6_to_UdJH(F0, F2, F4, F6):
    """
    Given :math:`F_0`, :math:`F_2`, :math:`F_4` and :math:`F_6`,
    return :math:`U_d` and :math:`J_H`, for :math:`f`-orbitals.

    Parameters
    ----------
    F0: float
        Slater integral :math:`F_0`.
    F2: float
        Slater integral :math:`F_2`.
    F4: float
        Slater integral :math:`F_4`.
    F6: float
        Slater integral :math:`F_6`.

    Returns
    -------
    Ud: float
        Coulomb interaction :math:`U_d`.
    JH: float
        Hund's coupling :math:`J_H`.
    """

    Ud = F0
    JH = (286 * F2 + 195 * F4 + 250 * F6) / 6435.0
    return Ud, JH


def CT_imp_bath(U_dd, Delta, n):
    """
    Compute energies of the impurity and bath for an
    Anderson impurity or charge-transfer model
    appropriate for a :math:`d`-shell transition metal compound.

    Parameters
    ----------
    U_dd: float
        Coulomb interaction :math:`U_{dd}`
    Delta: float
        Charge-transfer energy :math:`\\Delta`
    n : integer
        Number of electrons in the :math:`d`-shell

    Returns
    -------
    E_d : float
        Energy of the impurity states :math:`E_d`
    E_L : float
        Energy of the bath states :math:`E_L`

    Notes
    -----
    We credit our approach to Maurits Hakverkort,
    Heidelberg University.
    We define the state with a full set of bath orbitals to be zero
    energy and write the energy levels using the same definitions
    as [1]_ [2]_ [3]_.

    * :math:`d^{n}L^{10}` has energy :math:`0`

    * :math:`d^{n+1}L^9` has energy :math:`\\Delta`

    * :math:`d^{n+2}L^8` has energy :math:`2\\Delta + U_{dd}`

    Using this we can write and solve three linear equations to get
    :math:`E_d` and :math:`E_L` the energies of the impurity and bath.

       .. math::
           \\begin{aligned}
           10 E_L + n     E_d + n(n-1) \\frac{U_{dd}}{2} &= 0 \\\\
            9 E_L + (n+1) E_d + (n+1)n \\frac{U_{dd}}{2} &= \\Delta \\\\
            8 E_L + (n+2) E_d + (n+1)(n+2) \\frac{U_{dd}}{2}
            &= 2\\Delta + U_{dd}
           \\end{aligned}

    The solutions are:

       .. math::
           \\begin{aligned}
           E_d &=
           \\frac{10 \\Delta - n (19 +  n) U_{dd}/2}{10 + n}  \\\\
           E_L &=
           \\frac{n ((1+n) U_{dd}/2-\\Delta)}{10 + n}
           \\end{aligned}.

    References
    ----------
    .. [1] J. Zaanen, G. A. Sawatzky, and J. W. Allen
           `Phys. Rev. Lett. 55, 418 (1985) <https://doi.org/10.1103/PhysRevLett.55.418>`_
    .. [2] Maurits Haverkort et al.,
           `Phys. Rev. B 85, 165113 (2012) <https://doi.org/10.1103/PhysRevB.85.165113>`_
    .. [3] A. E. Bocquet et al.,
           `Phys. Rev. B 53, 1161 (1996) <https://doi.org/10.1103/PhysRevB.53.1161>`_
    """
    E_d = (10*Delta - n*(19 + n)*U_dd/2)/(10 + n)
    E_L = n*((1 + n)*U_dd/2-Delta)/(10 + n)
    return E_d, E_L


def CT_imp_bath_core_hole(U_dd, U_pd, Delta, n):
    """
    Compute energies of the impurity and bath for an
    Anderson impurity or charge-transfer model
    appropriate for a :math:`d`-shell transition metal compound
    with a core hole.

    Parameters
    ----------
    U_dd: float
        Coulomb interaction :math:`U_{dd}`
    U_pd: float
        Coulomb interaction :math:`U_{pd}`
    Delta: float
        Charge-transfer energy :math:`\\Delta`
    n : integer
        Number of electrons in the :math:`d`-shell :math:`n`

    Returns
    -------
    E_dc : float
        Energy of the impurity states :math:`E_{dc}` with
        a core hole
    E_Lc : float
        Energy of the bath states :math:`E_{Lc}` with a core
        hole
    E_p : float
        Energy of the core hole state :math:`E_\\textrm{E_p}`

    Notes
    -----
    We credit our approach to Maurits Hakverkort,
    Heidelberg University.
    We define the state with a full set of bath orbitals to be zero
    energy and write the energy levels using the same definitions as
    [1]_ [2]_ [3]_.

    * :math:`2p^5 d^{n}L^{10}` has energy :math:`0`

    * :math:`2p^5 d^{n+1}L^9` has energy :math:`\\Delta + U_{dd} - U_{pd}`

    * :math:`2p^5 d^{n+2}L^8` has energy :math:`2\\Delta + 3 U_{dd} - 2 U_{pd}`

    Using this we can write and solve linear equations to get
    :math:`E_{dc}`, :math:`E_{Lc}` and :math:`E_p` the energies of the
    impurity and bath with a core hole and the energy of the core hole.

       .. math::
           \\begin{aligned}
               6 E_p + 10 E_{Lc} +  n    E_{dc} + n(n-1) \\frac{U_{dd}}{2}
               + 6 n     U_{pd} &= 0    \\\\
               6 E_p +  9 E_{Lc} + (n+1) E_{dc} + (n+1)n \\frac{U_{dd}}{2}
               + 6 (n+1) U_{pd} &= \\Delta  \\\\
               6 E_p +  8 E_{Lc} + (n+2) E_{d} + (n+1)(n+2) \\frac{U_{dd}}{2}
               + 6 (n+2) U_{pd} &= 2 \\Delta+U_{dd} \\\\
               5 E_p + 10 E_{Lc} + (n+1) E_{d} + (n+1) n \\frac{U_{dd}}{2}
               + 5 (n+1) U_{pd} &= 0 \\\\
               5 E_p +  9 E_{Lc} + (n+2) E_{dc} + (n+2)(n+1) \\frac{U_{dd}}{2}
               + 5 (n+2) U_{pd} &= \\Delta+U_{dd}-U_{pd} \\\\
               5 E_p +  8 E_{Lc} + (n+3) E_{dc} + (n+3)(n+2) \\frac{U_{dd}}{2}
               + 5 (n+3) U_{pd} &= 2 \\Delta+3 U_{dd}-2 U_{pd}
           \\end{aligned}

    The solutions are:

       .. math::
           \\begin{aligned}
              E_{dc} &= \\frac{10 \\Delta - n (31+n) \\frac{U_{dd}}{2}-90 U_{pd}}{16+n} \\\\
              E_{Lc} &= \\frac{(1+n) (n \\frac{U_{dd}}{2}+6*U_{pd})-(6+n) \\Delta}{16+n} \\\\
              E_p    &= \\frac{10 \\Delta + (1+n)(n\\frac{U_{dd}}{2}-(10+n)*U_{pd}}{16+n}
           \\end{aligned}

    References
    ----------
    .. [1] J. Zaanen, G. A. Sawatzky, and J. W. Allen
           `Phys. Rev. Lett. 55, 418 (1985) <https://doi.org/10.1103/PhysRevLett.55.418>`_
    .. [2] Maurits Haverkort et al.,
           `Phys. Rev. B 85, 165113 (2012) <https://doi.org/10.1103/PhysRevB.85.165113>`_
    .. [3] A. E. Bocquet et al.,
           `Phys. Rev. B 53, 1161 (1996) <https://doi.org/10.1103/PhysRevB.53.1161>`_
    """
    E_dc = (10*Delta - n*(31 + n)*U_dd/2 - 90*U_pd) / (16 + n)
    E_Lc = ((1 + n)*(n*U_dd/2 + 6*U_pd) - (6 + n)*Delta) / (16 + n)
    E_p = (10*Delta + (1 + n)*(n*U_dd/2 - (10 + n)*U_pd)) / (16 + n)
    return E_dc, E_Lc, E_p


def info_atomic_shell():
    """
    Return a dict to describe the information of atomic shell.
    The key is a string of the shell's name, the value is a 2-elments tuple,
    the first value is the orbital angular moment number and
    the second value is the dimension of the Hilbert space.

    Returns
    -------
    info: dict
        The dict describing info of atomic shell.
    """

    info = {'s':   (0, 2),
            'p':   (1, 6),
            'p12': (1, 2),
            'p32': (1, 4),
            't2g': (2, 6),
            'd':   (2, 10),
            'd32': (2, 4),
            'd52': (2, 6),
            'f':   (3, 14),
            'f52': (3, 6),
            'f72': (3, 8)
            }

    return info


def case_to_shell_name(case):
    """
    Return the shell names for different cases.

    Parameters
    ----------
    case: string
        A string describing the shells included.

    Returns
    -------
    shell_name: tuple of one or two strings
        The name of shells.

    Examples
    --------
    >>> import edrixs
    >>> edrixs.case_to_shell_name('d')
    ('d',)
    >>> edrixs.case_to_shell_name('t2gp32')
    ('t2g', 'p32')
    """

    shell = ['s', 'p', 'p12', 'p32', 't2g', 'd', 'd32', 'd52', 'f', 'f52', 'f72']

    shell_name = {}
    for str1 in shell:
        shell_name[str1] = (str1,)

    for str1 in shell:
        for str2 in shell:
            shell_name[str1+str2] = (str1, str2)

    return shell_name[case.strip()]


def edge_to_shell_name(edge_name, with_main_qn=False):
    """
    Given edge name, return shell name.
    If one needs to include both spin-orbit split edges of one shell,
    one can use string, for example, 'L23' means both L2 and L3 edges
    are considered in the calculations, and the shell name will be 'p'.

    Parameters
    ----------
    edge_name: string
        Standard edge name.
    with_main_qn: logical
        If true, the shell name will include the main quantum number.
    Returns
    -------
    shell_name: string
        Shell name.

        - with_main_qn=True, shell_name will include the main quantum number,
          for example, *2p*

        - with_main_qn=False, shell_name will not include the main quantum number,
          for example, *p*
    """
    shell_name = {
        'K': ('s', '1s'),
        'L1': ('s', '2s'),
        'L2': ('p12', '2p12'),
        'L3': ('p32', '2p32'),
        'L23': ('p', '2p'),
        'M1': ('s', '3s'),
        'M2': ('p12', '3p12'),
        'M3': ('p32', '3p32'),
        'M23': ('p', '3p'),
        'M4': ('d32', '3d32'),
        'M5': ('d52', '3d52'),
        'M45': ('d', '3d'),
        'N1': ('s', '4s'),
        'N2': ('p12', '4p12'),
        'N3': ('p32', '4p32'),
        'N23': ('p', '4p'),
        'N4': ('d32', '4d32'),
        'N5': ('d52', '4d52'),
        'N45': ('d', '4d'),
        'N6': ('f52', '4f52'),
        'N7': ('f72', '4f72'),
        'N67': ('f', '4f'),
        'O1': ('s', '5s'),
        'O2': ('p12', '5p12'),
        'O3': ('p32', '5p32'),
        'O23': ('p', '5p'),
        'O4': ('d32', '5d32'),
        'O5': ('d52', '5d52'),
        'O45': ('d', '5d'),
        'P1': ('s', '6s'),
        'P2': ('p12', '6p12'),
        'P3': ('p32', '6p32'),
        'P23': ('p', '6p')
    }

    if with_main_qn:
        return shell_name[edge_name.strip()][1]
    else:
        return shell_name[edge_name.strip()][0]


def slater_integrals_name(shell_name, label=None):
    """
    Given shell names, return the names of required Slater integrals.
    The order of these names in the return list is according to the convention:

    - For 1 shell: [FX_11]

    - For 2 shells: [FX_11, FX_12, GX_12, FX_22]

    - For 3 shells: [FX_11, FX_12, GX_12, FX_22, FX_13, GX_13, FX_23, GX_23, FX_33]

    where, X=0, 2, 4, 6, or X=1, 3, 5, X shoule be in ascending order.

    Parameters
    ----------
    shell_name: tuple of strings
        Name of shells. Its length should be less and equal than 3.
    label: tuple of strings
        Label of shells, same shape as shell_name.

        If not provided, label will be set to

        - label=(1,) for one shell, or

        - label=(1,2) for two shells, or

        - label=(1,2,3) for three shells

    Returns
    -------
    res: list of strings
        Names of Slater integrals.
    """
    info = info_atomic_shell()
    # one shell
    if len(shell_name) == 1:
        res = []
        l1 = info[shell_name[0]][0]
        if label is not None:
            x = label[0]
        else:
            x = '1'
        res.extend(['F' + str(i) + '_' + x + x for i in range(0, 2*l1+1, 2)])
    elif len(shell_name) == 2:
        res = []
        l1 = info[shell_name[0]][0]
        l2 = info[shell_name[1]][0]
        if label is not None:
            x, y = label[0], label[1]
        else:
            x, y = '1', '2'
        res.extend(['F' + str(i) + '_' + x + x for i in range(0, 2*l1+1, 2)])
        res.extend(['F' + str(i) + '_' + x + y for i in range(0, min(2*l1, 2*l2)+1, 2)])
        res.extend(['G' + str(i) + '_' + x + y for i in range(abs(l1-l2), l1+l2+1, 2)])
        res.extend(['F' + str(i) + '_' + y + y for i in range(0, 2*l2+1, 2)])
    elif len(shell_name) == 3:
        res = []
        l1 = info[shell_name[0]][0]
        l2 = info[shell_name[1]][0]
        l3 = info[shell_name[2]][0]
        if label is not None:
            x, y, z = label[0], label[1], label[2]
        else:
            x, y, z = '1', '2', '3'
        res.extend(['F' + str(i) + '_' + x + x for i in range(0, 2*l1+1, 2)])
        res.extend(['F' + str(i) + '_' + x + y for i in range(0, min(2*l1, 2*l2)+1, 2)])
        res.extend(['G' + str(i) + '_' + x + y for i in range(abs(l1-l2), l1+l2+1, 2)])
        res.extend(['F' + str(i) + '_' + y + y for i in range(0, 2*l2+1, 2)])
        res.extend(['F' + str(i) + '_' + x + z for i in range(0, min(2*l1, 2*l3)+1, 2)])
        res.extend(['G' + str(i) + '_' + x + z for i in range(abs(l1-l3), l1+l3+1, 2)])
        res.extend(['F' + str(i) + '_' + y + z for i in range(0, min(2*l2, 2*l3)+1, 2)])
        res.extend(['G' + str(i) + '_' + y + z for i in range(abs(l2-l3), l2+l3+1, 2)])
        res.extend(['F' + str(i) + '_' + z + z for i in range(0, 2*l3+1, 2)])
    else:
        raise Exception("Not implemented for this case: ", shell_name)

    return res


def get_atom_data(atom, v_name, v_noccu, edge=None, trans_to_which=1, label=None):
    """
    Return Slater integrals, spin-orbit coupling strength, edge energy and core hole
    life-time broadening for an atom with given type of valence shells, occupancy of
    valence shells and the x-ray resonant edge.

    The Slater integrals and spin-orbit coupling are calculated by Cowan's code
    (https://github.com/mretegan/atomic-parameters) based on Hartree-Fock
    approximation. These numbers are just initial guess and serve as a start point,
    you need to rescale them to reproduce your XAS or RIXS spectra.

    NOTE: F0 is not calculated by Cowan's code, they are all set to be zero, you need to
    set F0 by yourself.

    Parameters
    ----------
    atom: string
        Name of atom, for example, 'Ni', 'Cu', 'Ir', 'U'.
    v_name: a string or a tuple of one or two strings
        Names of one or two valence shells. Set it to a single string if there is only
        one valence shell, or set it to a tuple of two strings if there are two valence
        shells. For example, *v_name='3d'*, *v_name=('3d', '4p')*
    v_noccu: a int or a tuple of one or two ints
        Occupancy of valence shells before x-ray absorption (without a core hole).
        Set it to a single int if there is only one valence shell, or set it to a
        tuple of two ints if there are two valence shells. The order should be
        consistent with that in v_name.
        For example, *v_name='3d', v_noccu=8*, *v_name=('3d', '4p'), v_noccu=(8, 0)*.
    edge: a string
        X-ray resonant edge. If edge is not None, both the information of
        the initial (without a core hole) and intermediate (with a core hole) Hamiltonians
        will be returned, otherwise, only the information of the initial Hamiltonian
        will be returned. It can be,

        - K, L1, L2, L3, M1, M2, M3, M4, M5, N1, N2, N3, N4, N5,
          N6, N7, O1, O2, O3, O4, O5, P1, P2, P3

        - L23 (both L2 and L3), M23 (both M2 and M3), M45 (both M4 and M5)
          N23 (both N2 and N3), N45 (both N4 and N5), N67 (both N6 and N7)
          O23 (both O2 and O3), O45 (both O4 and O5), P23 (both P2 and P3)
    trans_to_which: int
        If there are two valence shells, this variable is used to indicate which
        valence shell the photon transition happens.

        - 1: to the first valence shell

        - 2: to the second valence shell
    label: a string or tuple of strings
        User-defined symbols to label the names of Slater integrals.

        - one element, for the valence shell

        - two elements, 1st (2nd)-element for the 1st (2nd) valence shells or
          1st-element for the 1st valence shell and the 2nd one for the core shell.

        - three elements, the 1st (2nd)-element is for the 1st (2nd) valence shell,
          the 3rd one is for the core shell


    Returns
    -------
    res:  dict
        Atomic data. A dict likes this if edge is None:

        .. code-block:: python

           res={
               'slater_i': [(name of slater integrals, vaule of slater integrals), (), ...],
               'v_soc': [list of SOC for each atomic shell in v_name],
           }

        otherwise,

        .. code-block:: python

           res={
               'slater_i': [(name of slater integrals, value of slater integrals), (), ...],
               'slater_n': [(name of slater integrals, value of slater integrals), (), ...],
               'v_soc': [list of SOC for each atomic shell in v_name],
               'c_soc': SOC of core shell,
               'edge_ene': [list of edge energy, 2 elements if edge is any of
               "L23, M23, M45, N23, N45, N67, O23, O45, P23", othewise only 1 element],
               'gamma_c': [list of core hole life-time broadening, 2 elements if edge is any of
               "L23, M23, M45, N23, N45, N67, O23, O45, P23", othewise only 1 element],
           }


    Examples
    --------
    >>> import edrixs
    >>> res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8)
    >>> res
    {'slater_i': [('F0_11', 0.0), ('F2_11', 13.886), ('F4_11', 8.67)],
     'v_soc_i': [0.112]}
    >>> name, slater = [list(i) for i in zip(*res['slater_i'])]
    >>> name
    ['F0_11', 'F2_11', 'F4_11']
    >>> slater
    [0.0, 13.886, 8.67]
    >>> slater = [i * 0.8 for i in slater]
    >>> slater[0] = edrixs.get_F0('d', slater[1], slater[2])
    >>> slater
    [0.5728507936507936, 11.1088, 6.936]

    >>> import edrixs
    >>> res=edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8, edge='L3', label=('d', 'p'))
    >>> res
    {'slater_i': [('F0_dd', 0.0), ('F2_dd', 12.234), ('F4_dd', 7.598)],
     'v_soc_i': [0.083],
     'slater_n': [('F0_dd', 0.0), ('F2_dd', 12.234), ('F4_dd', 7.598),
      ('F0_dp', 0.0), ('F2_dp', 7.721), ('G1_dp', 5.787), ('G3_dp', 3.291),
      ('F0_pp', 0.0), ('F2_pp', 0.0)],
     'v_soc_n': [0.102],
     'c_soc': 11.507,
     'edge_ene': [852.7],
     'gamma_c': [0.275]}
    >>> name_i, slater_i = [list(i) for i in zip(*res['slater_i'])]
    >>> name_n, slater_n = [list(i) for i in zip(*res['slater_n'])]
    >>> name_i
    ['F0_dd', 'F2_dd', 'F4_dd']
    >>> slater_i
    [0.0, 12.234, 7.598]
    >>> name_n
    ['F0_dd', 'F2_dd', 'F4_dd', 'F0_dp', 'F2_dp', 'G1_dp', 'G3_dp', 'F0_pp', 'F2_pp']
    >>> slater_n
    [0.0, 12.234, 7.598, 0.0, 7.721, 5.787, 3.291, 0.0, 0.0]
    >>> slater_n[0] = edrixs.get_F0('d', slater_n[1], slater_n[2])
    >>> slater_n[3] = edrixs.get_F0('dp', slater_n[5], slater_n[6])
    >>> slater_n
    [0.6295873015873016, 12.234, 7.598, 0.5268428571428572, 7.721, 5.787, 3.291, 0.0, 0.0]

    >>> import edrixs
    >>> import collections
    >>> res=edrixs.get_atom_data('Ni', v_name=('3d', '4p'), v_noccu=(8, 0),
    ...     edge='K', trans_to_which=2, label=('d', 'p', 's'))
    >>> res
    {'slater_i': [('F0_dd', 0.0), ('F2_dd', 12.234), ('F4_dd', 7.598),
      ('F0_dp', 0.0), ('F2_dp', 0.0), ('G1_dp', 0.0), ('G3_dp', 0.0),
      ('F0_pp', 0.0), ('F2_pp', 0.0)],
     'v_soc_i': [0.083, 0.0],
     'slater_n': [('F0_dd', 0.0), ('F2_dd', 13.851), ('F4_dd', 8.643),
      ('F0_dp', 0.0), ('F2_dp', 2.299), ('G1_dp', 0.828), ('G3_dp', 0.713),
      ('F0_pp', 0.0), ('F2_pp', 0.0),
      ('F0_ds', 0.0), ('G2_ds', 0.079),
      ('F0_ps', 0.0), ('G1_ps', 0.194),
      ('F0_ss', 0.0)],
     'v_soc_n': [0.113, 0.093],
     'c_soc': 0.0,
     'edge_ene': [8333.0],
     'gamma_c': [0.81]}
    >>> slat_i = collections.OrderedDict(res['slater_i'])
    >>> slat_n = collections.OrderedDict(res['slater_n'])
    >>> list(slat_i.keys())
    ['F0_dd', 'F2_dd', 'F4_dd', 'F0_dp', 'F2_dp', 'G1_dp', 'G3_dp', 'F0_pp', 'F2_pp']
    >>> list(slat_i.values())
    [0.0, 12.234, 7.598, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    >>> list(slat_n.keys())
    ['F0_dd', 'F2_dd', 'F4_dd', 'F0_dp', 'F2_dp', 'G1_dp', 'G3_dp', 'F0_pp', 'F2_pp',
     'F0_ds', 'G2_ds', 'F0_ps', 'G1_ps', 'F0_ss']
    >>> list(slat_n.values())
    [0.0, 13.851, 8.643, 0.0, 2.299, 0.828, 0.713, 0.0, 0.0,
     0.0, 0.079, 0.0, 0.194, 0.0]
    """
    c_norb = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
    atom = atom.strip()
    avail_atoms = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                   'Re', 'Os', 'Ir',
                   'Sm',
                   'U', 'Pu']
    avail_shells = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f',
                    '5s', '5p', '5d', '5f', '6s', '6p', '6d']

    if atom not in avail_atoms:
        raise Exception("Atom data is Not available for this atom: ", atom)

    if not isinstance(v_name, (list, tuple)):
        v_name = (v_name,)
    if not isinstance(v_noccu, (list, tuple)):
        v_noccu = (v_noccu,)

    if label is not None:
        if not isinstance(label, (list, tuple)):
            label = (label,)

    if len(v_name) != len(v_noccu):
        raise Exception("The shape of v_name is not same as noccu")

    for ishell in v_name:
        if ishell not in avail_shells:
            raise Exception("Not available for this shell: ", ishell)

    data_pkg = 'edrixs.atom_data'
    fname = f'{atom}.json'
    with resources.open_text(data_pkg, fname) as f:
        atom_dict = json.load(f)

    res = {}
    shell_name = []
    for name in v_name:
        shell_name.append(name[-1])
    if label is not None:
        my_label = label[0:len(shell_name)]
    else:
        my_label = None
    slater_name = slater_integrals_name(shell_name, label=my_label)

    if len(v_name) == 1:
        case = v_name[0] + str(v_noccu[0])
    else:
        case = v_name[0] + str(v_noccu[0]) + '_' + v_name[1] + str(v_noccu[1])
    if case not in atom_dict:
        raise Exception("This configuration is not available in atom_data", case)

    nslat = len(slater_name)
    slater_i = [0.0] * nslat
    tmp = atom_dict[case]['slater']
    slater_i[0:len(tmp)] = tmp
    res['slater_i'] = list(zip(slater_name, slater_i))
    res['v_soc_i'] = atom_dict[case]['soc']

    if edge is not None:
        edge = edge.strip()
        edge_name = edge_to_shell_name(edge, with_main_qn=True)
        shell_name = []
        for name in v_name:
            shell_name.append(name[-1])
        shell_name.append(edge_name[1:2])
        if label is not None:
            my_label = label[0:len(shell_name)]
        else:
            my_label = None
        slater_name = slater_integrals_name(shell_name, label=my_label)

        if len(v_name) == 1:
            case = (v_name[0] + str(v_noccu[0]+1) + '_' +
                    edge_name[0:2] + str(c_norb[edge_name[1]]-1))
        else:
            if trans_to_which == 1:
                case = (v_name[0] + str(v_noccu[0]+1) + '_' +
                        v_name[1] + str(v_noccu[1]) + '_' +
                        edge_name[0:2] + str(c_norb[edge_name[1]]-1))
            else:
                case = (v_name[0] + str(v_noccu[0]) + '_' +
                        v_name[1] + str(v_noccu[1]+1) + '_' +
                        edge_name[0:2] + str(c_norb[edge_name[1]]-1))
        if case not in atom_dict:
            raise Exception("This configuration is currently not available in atom_data", case)

        nslat = len(slater_name)
        slater_n = [0.0] * nslat
        tmp = atom_dict[case]['slater']
        slater_n[0:len(tmp)] = tmp
        res['slater_n'] = list(zip(slater_name, slater_n))

        res['v_soc_n'] = atom_dict[case]['soc'][0:-1]
        res['c_soc'] = atom_dict[case]['soc'][-1]

        res['edge_ene'] = atom_dict[edge]['ene']
        res['gamma_c'] = [i / 2 for i in atom_dict[edge]['gamma']]

    return res


def rescale(old_list, scale=None):
    """
    Rescale a 1d list.

    Parameters
    ----------
    old_list: 1d list of numerical values
        Input list.
    scale: a list of tuples
        The scaling factors.

    Returns
    -------
    new_list: 1d list
        The rescaled list.
    """

    new_list = [i for i in old_list]
    if scale is not None:
        for pos, val in zip(scale[0], scale[1]):
            new_list[pos] = new_list[pos] * val
    return new_list
