__all__ = ['get_ladd', 'get_lminus', 'get_lx', 'get_ly', 'get_lz', 'get_orb_momentum',
           'get_pauli', 'get_sx', 'get_sy', 'get_sz', 'get_spin_momentum', 'euler_to_rmat',
           'rmat_to_euler', 'where_is_angle', 'dmat_spinor', 'zx_to_rmat', 'get_wigner_dmat',
           'cf_cubic_d', 'cf_tetragonal_d', 'cf_square_planar_d', 'cf_trigonal_t2g',
           'cf_trigonal_d']

import numpy as np
from .basis_transform import cb_op, tmat_r2c, tmat_trig2r


def get_ladd(ll, ispin=False):
    """
    Get the matrix form of the raising operator :math:`l^+` in the
    complex spherical harmonics basis

    .. math::

        l^+|ll,m> = \\sqrt{(ll-m)(ll+m+1)} |ll,m+1>

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    ladd: 2d complex array
        The matrix form of :math:`l^+`.

        If ispin=True, the dimension will be :math:`2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`(2ll+1) \\times (2ll+1)`.
    """

    norbs = 2 * ll + 1
    ladd = np.zeros((norbs, norbs), dtype=np.complex128)
    cone = np.complex128(1.0 + 0.0j)
    for m in range(-ll, ll):
        ladd[ll + m + 1, ll + m] = np.sqrt((ll - m) * (ll + m + 1.0) * cone)
    if ispin:
        ladd_spin = np.zeros((2 * norbs, 2 * norbs), dtype=np.complex128)
        ladd_spin[0:2 * norbs:2, 0:2 * norbs:2] = ladd
        ladd_spin[1:2 * norbs:2, 1:2 * norbs:2] = ladd
        return ladd_spin
    else:
        return ladd


def get_lminus(ll, ispin=False):
    """
    Get the matrix form of the lowering operator :math:`l^-` in the
    complex spherical harmonics basis

    .. math::

        l^-|ll,m> = \\sqrt{(ll+m)(ll-m+1)} |ll,m-1>

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    lminus: 2d complex array
        The matrix form of :math:`l^-`.

        If ispin=True, the dimension will be :math:`2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`(2ll+1) \\times (2ll+1)`.
    """

    norbs = 2 * ll + 1
    lminus = np.zeros((norbs, norbs), dtype=np.complex128)
    cone = np.complex128(1.0 + 0.0j)
    for m in range(-ll + 1, ll + 1):
        lminus[ll + m - 1, ll + m] = np.sqrt((ll + m) * (ll - m + 1.0) * cone)
    if ispin:
        lminus_spin = np.zeros((2 * norbs, 2 * norbs), dtype=np.complex128)
        lminus_spin[0:2 * norbs:2, 0:2 * norbs:2] = lminus
        lminus_spin[1:2 * norbs:2, 1:2 * norbs:2] = lminus
        return lminus_spin
    else:
        return lminus


def get_lx(ll, ispin=False):
    """
    Get the matrix form of the orbital angular momentum
    operator :math:`l_x` in the complex spherical harmonics basis,

    .. math::

        l_x = \\frac{1}{2} (l^+ + l^-)

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    lx: 2d complex array
        The matrix form of :math:`l_x`.

        If ispin=True, the dimension will be :math:`2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`(2ll+1) \\times (2ll+1)`.
    """

    if ispin:
        return (get_ladd(ll, True) + get_lminus(ll, True)) / 2.0
    else:
        return (get_ladd(ll) + get_lminus(ll)) / 2.0


def get_ly(ll, ispin=False):
    """
    Get the matrix form of the orbital angular momentum
    operator :math:`l_y` in the complex spherical harmonics basis,

    .. math::

        l_y = \\frac{-i}{2} (l^+ - l^-)

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    ly: 2d complex array
        The matrix form of :math:`l_y`.

        If ispin=True, the dimension will be :math:`2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`(2ll+1) \\times (2ll+1)`.
    """

    if ispin:
        return (get_ladd(ll, True) - get_lminus(ll, True)) / 2.0 * -1j
    else:
        return (get_ladd(ll) - get_lminus(ll)) / 2.0 * -1j


def get_lz(ll, ispin=False):
    """
    Get the matrix form of the orbital angular momentum
    operator :math:`l_z` in the complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    lz: 2d complex array
        The matrix form of :math:`l_z`.

        If ispin=True, the dimension will be :math:`2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`(2ll+1) \\times (2ll+1)`.
    """
    norbs = 2 * ll + 1
    lz = np.zeros((norbs, norbs), dtype=np.complex128)
    for m in range(-ll, ll + 1):
        lz[ll + m, ll + m] = m
    if ispin:
        lz_spin = np.zeros((2 * norbs, 2 * norbs), dtype=np.complex128)
        lz_spin[0:2 * norbs:2, 0:2 * norbs:2] = lz
        lz_spin[1:2 * norbs:2, 1:2 * norbs:2] = lz
        return lz_spin
    else:
        return lz


def get_orb_momentum(ll, ispin=False):
    """
    Get the matrix form of the orbital angular momentum
    operator :math:`l_x, l_y, l_z` in the complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.
    ispin: logical
        Whether including spin or not (default: False).

    Returns
    -------
    res: 3d complex array
        The matrix form of

        - res[0]: :math:`l_x`

        - res[1]: :math:`l_y`

        - res[2]: :math:`l_z`

        If ispin=True, the dimension will be :math:`3 \\times 2(2ll+1) \\times 2(2ll+1)`,

        otherwise, it will be :math:`3 \\times (2ll+1) \\times (2ll+1)`.
    """
    norbs = 2 * ll + 1
    if ispin:
        res = np.zeros((3, 2*norbs, 2*norbs), dtype=np.complex128)
    else:
        res = np.zeros((3, norbs, norbs), dtype=np.complex128)
    res[0] = get_lx(ll, ispin)
    res[1] = get_ly(ll, ispin)
    res[2] = get_lz(ll, ispin)

    return res


def get_pauli():
    """
    Get the Pauli matrix

    Returns
    -------
    sigma: 3d complex array, shape=(3, 2, 2)

        - sigma[0] is :math:`\\sigma_x`,

        - sigma[1] is :math:`\\sigma_y`,

        - sigma[2] is :math:`\\sigma_z`,
    """

    sigma = np.zeros((3, 2, 2), dtype=np.complex128)
    sigma[0, 0, 1] = 1.0
    sigma[0, 1, 0] = 1.0
    sigma[1, 0, 1] = -1.0j
    sigma[1, 1, 0] = 1.0j
    sigma[2, 0, 0] = 1.0
    sigma[2, 1, 1] = -1.0

    return sigma


def get_sx(ll):
    """
    Get the matrix form of spin angular momentum operator :math:`s_x` in the
    complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Quantum number of orbital angular momentum.

    Returns
    -------
    sx: 2d complex array.
        Matrix form of :math:`s_x`, the dimension
        is :math:`2(2ll+1) \\times 2(2ll+1)`,

        Orbital order is: \\|-ll,up\\>, \\|-ll,down\\>, ...,
        \\|+ll, up\\>, \\|+ll,down\\>.
    """

    norbs = 2 * (2 * ll + 1)
    sx = np.zeros((norbs, norbs), dtype=np.complex128)
    sigma = get_pauli()
    for i in range(2 * ll + 1):
        sx[2 * i:2 * i + 2, 2 * i:2 * i + 2] = sigma[0, :, :] / 2.0

    return sx


def get_sy(ll):
    """
    Get the matrix form of spin angular momentum operator :math:`s_y` in the
    complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Quantum number of orbital angular momentum.

    Returns
    -------
    sy: 2d complex array.
        Matrix form of :math:`s_y`, the dimension
        is :math:`2(2ll+1) \\times 2(2ll+1)`, spin order is:

        Orbital order is: \\|-ll,up\\>, \\|-ll,down\\>, ...,
        \\|+ll, up\\>, \\|+ll,down\\>
    """

    norbs = 2 * (2 * ll + 1)
    sy = np.zeros((norbs, norbs), dtype=np.complex128)
    sigma = get_pauli()
    for i in range(2 * ll + 1):
        sy[2 * i:2 * i + 2, 2 * i:2 * i + 2] = sigma[1, :, :] / 2.0

    return sy


def get_sz(ll):
    """
    Get the matrix form of spin angular momentum operator :math:`s_z` in the
    complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Quantum number of orbital angular momentum.

    Returns
    -------
    sz: 2d complex array.
        Matrix form of :math:`s_z`, the dimension
        is :math:`2(2ll+1) \\times 2(2ll+1)`.

        Orbital order is: \\|-ll,up\\>, \\|-ll,down\\>, ...,
        \\|+ll, up\\>, \\|+ll,down\\>
    """

    norbs = 2 * (2 * ll + 1)
    sz = np.zeros((norbs, norbs), dtype=np.complex128)
    sigma = get_pauli()
    for i in range(2 * ll + 1):
        sz[2 * i:2 * i + 2, 2 * i:2 * i + 2] = sigma[2, :, :] / 2.0

    return sz


def get_spin_momentum(ll):
    """
    Get the matrix form of the spin angular momentum
    operator :math:`s_x, s_y, s_z` in the complex spherical harmonics basis.

    Parameters
    ----------
    ll: int
        Orbital angular momentum number.

    Returns
    -------
    res: 3d complex array
        The matrix form of

        - res[0]: :math:`s_x`

        - res[1]: :math:`s_y`

        - res[2]: :math:`s_z`

        the dimension is :math:`3 \\times 2(2ll+1) \\times 2(2ll+1)`,

        Orbital order is: \\|-ll,up\\>, \\|-ll,down\\>, ...,
        \\|+ll, up\\>, \\|+ll,down\\>
    """
    norbs = 2 * (2 * ll + 1)
    res = np.zeros((3, norbs, norbs), dtype=np.complex128)
    res[0] = get_sx(ll)
    res[1] = get_sy(ll)
    res[2] = get_sz(ll)

    return res


def euler_to_rmat(alpha, beta, gamma):
    """
    Given Euler angle: :math:`\\alpha, \\beta, \\gamma`,
    generate the :math:`3 \\times 3` rotational matrix :math:`R`.

    Parameters
    ----------
    alpha: float
        Euler angle, in radian, [0, :math:`2\\pi`]
    beta: float
        Euler angle, in radian, [0, :math:`\\pi`]
    gamma: float
        Euler angle, in radian, [0, :math:`2\\pi`]

    Returns
    -------
    rmat: 2d float array
        The :math:`3 \\times 3` rotational matrix.
    """

    rmat = np.zeros((3, 3), dtype=np.float64)
    rmat[0, 0] = np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.sin(alpha) * np.sin(gamma)
    rmat[0, 1] = -np.sin(gamma) * np.cos(alpha) * np.cos(beta) - np.sin(alpha) * np.cos(gamma)
    rmat[0, 2] = np.cos(alpha) * np.sin(beta)
    rmat[1, 0] = np.sin(alpha) * np.cos(beta) * np.cos(gamma) + np.cos(alpha) * np.sin(gamma)
    rmat[1, 1] = -np.sin(gamma) * np.sin(alpha) * np.cos(beta) + np.cos(alpha) * np.cos(gamma)
    rmat[1, 2] = np.sin(alpha) * np.sin(beta)
    rmat[2, 0] = -np.cos(gamma) * np.sin(beta)
    rmat[2, 1] = np.sin(beta) * np.sin(gamma)
    rmat[2, 2] = np.cos(beta)
    return rmat


def rmat_to_euler(rmat):
    """
    Given the :math:`3 \\times 3` rotational matrix :math:`R`, return the Euler
    angles: :math:`\\alpha, \\beta, \\gamma`.

    Parameters
    ----------
    rmat:  2d float array
        The :math:`3 \\times 3` rotational matrix :math:`R`.

    Returns
    -------
    alpha: float
        Euler angle :math:`\\alpha` in radian, [0, :math:`2\\pi`].
    beta: float
        Euler angle :math:`\\beta` in radian, [0, :math:`\\pi`].
    gamma: float
        Euler angle :math:`\\gamma` in radian, [0, :math:`2\\pi`]

    """

    if np.abs(rmat[2, 2]) < 1.0:
        beta = np.arccos(rmat[2, 2])

        cos_gamma = -rmat[2, 0] / np.sin(beta)
        sin_gamma = rmat[2, 1] / np.sin(beta)
        gamma = where_is_angle(sin_gamma, cos_gamma)

        cos_alpha = rmat[0, 2] / np.sin(beta)
        sin_alpha = rmat[1, 2] / np.sin(beta)
        alpha = where_is_angle(sin_alpha, cos_alpha)
    else:
        if rmat[2, 2] > 0:
            beta = 0.0
        else:
            beta = np.pi
        gamma = 0.0
        alpha = np.arccos(rmat[1, 1])
    return alpha, beta, gamma


def where_is_angle(sina, cosa):
    """
    Given sine and cosine of an angle :math:`\\alpha`, return the
    angle :math:`\\alpha` range from [0, :math:`2\\pi`].

    Parameters
    ----------
    sina:  float
        :math:`\\sin(\\alpha)`.
    cosa:  float
        :math:`\\cos(\\alpha)`.

    Returns
    -------
    alpha: float
        The angle :math:`\\alpha` in radian [0, :math:`2\\pi`].
    """

    if cosa > 1.0:
        cosa = 1.0
    elif cosa < -1.0:
        cosa = -1.0
    alpha = np.arccos(cosa)
    if sina < 0.0:
        alpha = 2.0 * np.pi - alpha
    return alpha


def dmat_spinor(alpha, beta, gamma):
    """
    Given three Euler angle: :math:`\\alpha, \\beta, \\gamma`,
    return the transformation
    matrix for :math:`\\frac{1}{2}`-spinor.

    Parameters
    ----------
    alpha: float
        Euler angle :math:`\\alpha` in radian [0, :math:`2\\pi`].
    beta: float
        Euler angle :math:`\\beta` in radian [0, :math:`\\pi`].
    gamma: float
        Euler angle :math:`\\gamma` in radian [0, :math:`2\\pi`].

    Returns
    -------
    dmat:  2d complex array
        The :math:`2 \\times 2` transformation matrix.
    """

    dmat = np.zeros((2, 2), dtype=np.complex128)
    dmat[0, 0] = np.exp(-(alpha + gamma) / 2.0 * 1j) * np.cos(beta / 2.0)
    dmat[0, 1] = -np.exp(-(alpha - gamma) / 2.0 * 1j) * np.sin(beta / 2.0)
    dmat[1, 0] = np.exp((alpha - gamma) / 2.0 * 1j) * np.sin(beta / 2.0)
    dmat[1, 1] = np.exp((alpha + gamma) / 2.0 * 1j) * np.cos(beta / 2.0)
    return dmat


def zx_to_rmat(z, x):
    """
    Given :math:`z` vector and :math:`x` vector, calculate :math:`y` vector
    which satisfies the right-hand Cartesian coordinate and normalize them to
    unit if needed, and then return the :math:`3 \\times 3`
    rotational matrix :math:`R`.

    Parameters
    ----------
    z: 1d float array
        The :math:`z` vector.
    x: 1d float array
        The :math:`x` vector.

    Returns
    -------
    rmat: 2d float array
        The :math:`3 \\times 3` rotational matrix :math:`R`.
    """

    z = np.array(z, dtype=np.float64)
    x = np.array(x, dtype=np.float64)
    xx = x / np.sqrt(np.dot(x, x))
    zz = z / np.sqrt(np.dot(z, z))
    yy = np.cross(zz, xx)

    rmat = np.zeros((3, 3), dtype=np.float64)
    rmat[:, 0] = xx
    rmat[:, 1] = yy
    rmat[:, 2] = zz

    return rmat


def get_wigner_dmat(quant_2j, alpha, beta, gamma):
    """
    Given quantum number and Euler angles, return the Wigner-D matrix.

    Parameters
    ----------
    quant_2j: int
        Twice of the quantum number j: 2j, for example, quant_2j=1 means j=1/2,
        quant_2j=2 means j=1
    alpha: float number
        The first Euler angle :math:`\\alpha` in radian [0, :math:`2\\pi`].
    beta: float number
        The second Euler angle :math:`\\beta` in radian [0, :math:`\\pi`].
    gamma: float number
        The third Euler angle :math:`\\gamma` in radian [0, :math:`2\\pi`].

    Returns
    -------
    result: 2d complex array, shape(quant_2j+1, quant_2j+1)
        The Wigner D-matrix.
        For :math:`j=1/2`, the orbital order is: +1/2 (spin up), -1/2 (spin down).
        For :math:`j>1/2`, the orbital order is: :math:`-j, -j+1, ..., +j`

    Examples
    --------
    >>> import edrixs
    spin-1/2 D-matrix
    >>> edrixs.get_wigner_dmat(1, 1, 2, 3)
    array([[-0.224845-0.491295j, -0.454649-0.708073j],
           [ 0.454649-0.708073j, -0.224845+0.491295j]])
    j=1 D-matrix
    >>> edrixs.get_wigner_dmat(2, 1, 2, 3)
    array([[-0.190816-0.220931j,  0.347398+0.541041j, -0.294663-0.643849j],
           [ 0.636536-0.090736j, -0.416147+0.j      , -0.636536-0.090736j],
           [-0.294663+0.643849j, -0.347398+0.541041j, -0.190816+0.220931j]])
    """

    from sympy.physics.quantum.spin import Rotation
    from sympy import N, S
    ndim = quant_2j + 1
    result = np.zeros((ndim, ndim), dtype=complex)
    # For j=1/2, we use different orbital order: first +1/2, then -1/2
    if quant_2j == 1:
        for i, mi in enumerate(range(quant_2j, -quant_2j-1, -2)):
            for j, mj in enumerate(range(quant_2j, -quant_2j-1, -2)):
                rot = Rotation.D(S(quant_2j)/2, S(mi)/2, S(mj)/2, alpha, beta, gamma)
                result[i, j] = N(rot.doit())
    # For j > 1/2, the order is -j, -j+1, ..., +j
    else:
        for i, mi in enumerate(range(-quant_2j, quant_2j+1, 2)):
            for j, mj in enumerate(range(-quant_2j, quant_2j+1, 2)):
                rot = Rotation.D(S(quant_2j)/2, S(mi)/2, S(mj)/2, alpha, beta, gamma)
                result[i, j] = N(rot.doit())

    return result


def cf_cubic_d(ten_dq):
    """
    Given 10Dq, return cubic crystal field matrix for d orbitals
    in the complex harmonics basis.

    Parameters
    ----------
    ten_dq: float scalar
        The splitting between :math:`eg` and :math:`t2g` orbitals.

    Returns
    -------
    cf: 2d complex array, shape=(10, 10)
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """

    tmp = np.zeros((5, 5), dtype=complex)
    tmp[0, 0] = 0.6 * ten_dq  # dz2
    tmp[1, 1] = -0.4 * ten_dq  # dzx
    tmp[2, 2] = -0.4 * ten_dq  # dzy
    tmp[3, 3] = 0.6 * ten_dq  # dx2-y2
    tmp[4, 4] = -0.4 * ten_dq  # dxy
    cf = np.zeros((10, 10), dtype=complex)
    cf[0:10:2, 0:10:2] = tmp
    cf[1:10:2, 1:10:2] = tmp

    cf[:, :] = cb_op(cf, tmat_r2c('d', True))

    return cf


def cf_tetragonal_d(ten_dq, d1, d3):
    """
    Given 10Dq, d1, d3, return tetragonal crystal field matrix for d orbitals
    in the complex harmonics basis.

    Parameters
    ----------
    ten_dq: float scalar
        Parameter used to label cubic crystal splitting.
    d1: float scalar
        Paramter used to label tetragonal splitting.
    d3: float scalar
        Paramter used to label tetragonal splitting.

    Returns
    -------
    cf: 2d complex array, shape=(10, 10)
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """

    dt = (3.0 * d3 - 4.0 * d1) / 35
    ds = (d3 + d1) / 7.0
    dq = ten_dq / 10.0

    tmp = np.zeros((5, 5), dtype=complex)
    tmp[0, 0] = 6 * dq - 2 * ds - 6 * dt  # d3z2-r2
    tmp[1, 1] = -4 * dq - 1 * ds + 4 * dt  # dzx
    tmp[2, 2] = -4 * dq - 1 * ds + 4 * dt  # dzy
    tmp[3, 3] = 6 * dq + 2 * ds - 1 * dt  # dx2-y2
    tmp[4, 4] = -4 * dq + 2 * ds - 1 * dt  # dxy

    cf = np.zeros((10, 10), dtype=complex)
    cf[0:10:2, 0:10:2] = tmp
    cf[1:10:2, 1:10:2] = tmp

    cf[:, :] = cb_op(cf, tmat_r2c('d', True))

    return cf


def cf_trigonal_t2g(delta):
    """
    Given delta, return trigonal crystal field matrix for t2g orbitals
    in the complex harmonics basis.

    Parameters
    ----------
    delta: float scalar
        Parameter used to label trigonal crystal splitting.

    Returns
    -------
    cf: 2d complex array, shape=(6, 6)
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """

    tmp = np.array([[0, delta, delta],
                    [delta, 0, delta],
                    [delta, delta, 0]])
    cf = np.zeros((6, 6), dtype=complex)
    cf[0:6:2, 0:6:2] += tmp
    cf[1:6:2, 1:6:2] += tmp
    cf[:, :] = cb_op(cf, tmat_r2c('t2g', True))

    return cf


def cf_trigonal_d(case, Delta1, Delta2):
    """
    Given Delta1 and Delta2, return trigonal crystal field matrix for d orbitals
    in the complex harmonics basis. See PHYSICAL REVIEW B 105, 245153 (2022) for the
    definition of the orbitals.

    Parameters
    ----------
    case: string
        Label for different systems.

        - '111': :math:`z` axis along diagonal of coordinate system
        - '100': :math:`z` axis vertical

    Delta1: float scalar
        The :math:`e^\\sigma_g` orbitals are Delta1 above the :math:`a^_{1g}`
    Delta1: float scalar
        The :math:`e^\\pi_g` orbitals are Delta2 below the :math:`a^_{1g}`

    Returns
    -------
    cf: 2d complex array
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """

    cf_trig = np.zeros((5, 5), dtype=complex)
    cf_trig[0, 0] = 0  # a1g
    cf_trig[1, 1] = -Delta2  # epi_g
    cf_trig[2, 2] = -Delta2  # epi_g
    cf_trig[3, 3] = Delta1  # esigma_g
    cf_trig[4, 4] = Delta1  # esigma_g

    cf_trig_spin = np.zeros((10, 10), dtype=complex)
    cf_trig_spin[::2, ::2] = cf_trig
    cf_trig_spin[1::2, 1::2] = cf_trig

    if case == '111':
        cf_real = cb_op(cf_trig_spin, tmat_trig2r('111', ispin=True))
        cf = cb_op(cf_real, tmat_r2c('d', ispin=True))
        return cf
    elif case == '100':
        cf_real = cb_op(cf_trig_spin, tmat_trig2r('100', ispin=True))
        cf = cb_op(cf_real, tmat_r2c('d', ispin=True))
        return cf
    else:
        raise Exception(f"case {case} not implemented."
                        "It should be either '111' or '110'")


def cf_square_planar_d(ten_dq, ds):
    """
    Given 10Dq, ds, return square planar crystal field matrix for d orbitals
    in the complex harmonics basis. This is the limit of strong tetragonal
    distortion with two axial ligands at infinity. Note that in this case the
    three parameters, ten_dq, ds, and dt, are no longer independent:
    dt = 2/35*ten_dq and the levels depend on only two parameters
    ten_dq and ds.

    Parameters
    ----------
    ten_dq: float scalar
        Parameter associated with eg-t2g splitting.
    ds: float scalar
        Paramter associated with splitting orbitals with
        z-components.

    Returns
    -------
    cf: 2d complex array, shape=(10, 10)
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """
    tmp = np.zeros((5, 5), dtype=complex)
    tmp[0, 0] = 9/35*ten_dq - 2*ds  # d3z2-r2
    tmp[1, 1] = -6/35*ten_dq - ds  # dzx
    tmp[2, 2] = -6/35*ten_dq - ds  # dzy
    tmp[3, 3] = 19/35*ten_dq + 2*ds  # dx2-y2
    tmp[4, 4] = -16/35*ten_dq + 2*ds  # dxy

    cf = np.zeros((10, 10), dtype=complex)
    cf[0:10:2, 0:10:2] = tmp
    cf[1:10:2, 1:10:2] = tmp

    cf[:, :] = cb_op(cf, tmat_r2c('d', True))

    return cf


def cf_tetragonal_t2g(ten_dq, d1, d3):
    """
    Given 10Dq, d1, d3, return tetragonal crystal field matrix for t2g orbitals
    in the complex harmonics basis.

    Parameters
    ----------
    ten_dq: float scalar
        Parameter used to label cubic crystal splitting.
    d1: float scalar
        Paramter used to label tetragonal splitting.
    d3: float scalar
        Paramter used to label tetragonal splitting.

    Returns
    -------
    cf: 2d complex array, shape=(6, 6)
        The matrix form of crystal field Hamiltonian in complex harmonics basis.
    """
    dt = (3.0 * d3 - 4.0 * d1) / 35
    ds = (d3 + d1) / 7.0
    dq = ten_dq / 10.0

    tmp = np.zeros((3, 3), dtype=complex)
    tmp[0, 0] = -4 * dq - 1 * ds + 4 * dt  # dxz
    tmp[1, 1] = -4 * dq - 1 * ds + 4 * dt  # dyz
    tmp[2, 2] = -4 * dq + 2 * ds - 1 * dt  # dxy

    cf = np.zeros((6, 6), dtype=complex)
    cf[0:6:2, 0:6:2] += tmp
    cf[1:6:2, 1:6:2] += tmp
    cf[:, :] = cb_op(cf, tmat_r2c('t2g', True))
    return cf
