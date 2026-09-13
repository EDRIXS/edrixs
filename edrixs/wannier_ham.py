__all__ = ['HR', 'KVec', 'SymKVec', 'UniKVec']

import numpy as np


class HR():
    """
    Class for post-process of the Wannier90 tight-binding (TB) Hamiltonian in real space.

    Parameters
    ----------
    nwann: int
        Number of Wannier orbitals.
    nrpt: int
        Number of :math:`r` points.
    irpt0: int
        Index of the point (0,0,0).
    rpts: float array
        The coordinates of :math:`r` points.
    deg_rpt: int array
        Degenerancy of :math:`r` points.
    hr: complex array
        Hamiltonian :math:`H(r)` from Wannier90.
    """

    def __init__(self, nwann, nrpt, irpt0, rpts, deg_rpt, hr):
        self.nwann = nwann
        self.nrpt = nrpt
        self.irpt0 = irpt0
        self.deg_rpt = deg_rpt
        self.rpts = rpts
        self.hr = hr

    @staticmethod
    def from_file(fname='wannier90_hr.dat'):
        """
        Generate a TB Hamiltonian from Wannier90 output file "case_hr.dat".

        Parameters
        ----------
        fname: str
            The file that contains the Wannier90 output file: "case_hr.dat".

        Returns
        -------
        HR: HR object
            A HR object.
        """

        with open(fname, 'r') as f:
            # skip the header (1 line)
            f.readline()
            nwann = int(f.readline().strip())
            nrpt = int(f.readline().strip())
            # the degeneracy of R points
            nline = nrpt // 15 + 1
            tmp = []
            for i in range(nline):
                tmp.extend(f.readline().strip().split())
            tmp = [int(item) for item in tmp]
            deg_rpt = np.array(tmp, dtype=int)
            # read hr for each r-point
            rpts = np.zeros((nrpt, 3), dtype=int)
            hr = np.zeros((nrpt, nwann, nwann), dtype=np.complex128)
            for i in range(nrpt):
                for j in range(nwann):
                    for k in range(nwann):
                        rx, ry, rz, hr_i, hr_j, hr_real, hr_imag = f.readline().strip().split()
                        rpts[i, :] = int(rx), int(ry), int(rz)
                        if int(rx) == 0 and int(ry) == 0 and int(rz) == 0:
                            irpt0 = i
                        hr[i, k, j] = np.float64(
                            hr_real) + np.float64(hr_imag) * 1j
            # construct the HR instance
            return HR(nwann, nrpt, irpt0, rpts, deg_rpt, hr)

    @staticmethod
    def copy_hr(other):
        """
        Copy instance of HR.

        Parameters
        ----------
        other: HR object
            A HR object to be copied.

        Returns
        -------
        HR: HR object
            Return a new HR object.
        """

        return HR(other.nwann,
                  other.nrpt,
                  other.irpt0,
                  np.copy(other.rpts),
                  np.copy(other.deg_rpt),
                  np.copy(other.hr))

    def get_hr0(self, ispin=False):
        """
        Return the on-site term :math:`H(r=0)`.

        Parameters
        ----------
        ispin: logical
            Whether to include spin degree of freedom or not (default: False).

        Returns
        -------
        hr: 2d complex array
            The on-site Hamiltonian.
        """

        if ispin:
            norbs = 2 * self.nwann
            hr0_spin = np.zeros((norbs, norbs), dtype=np.complex128)
            hr0_spin[0:norbs:2, 0:norbs:2] = self.hr[self.irpt0, :, :]
            hr0_spin[1:norbs:2, 1:norbs:2] = self.hr[self.irpt0, :, :]
            return hr0_spin
        else:
            return self.hr[self.irpt0, :, :]

    def get_hr(self, ispin):
        """
        Return the Hamiltonian of :math:`H(r)`.

        Parameters
        ----------
        ispin: logical
            Whether to include spin degree of freedom or not (default: False)

        Returns
        -------
        hr: 3d complex array
            The Hamiltonian :math:`H(r)`.
        """

        # with spin, spin order: up dn up dn up dn
        if ispin == 1:
            norbs = 2 * self.nwann
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:norbs:2, 0:norbs:2] = self.hr
            hr_spin[:, 1:norbs:2, 1:norbs:2] = self.hr
            return hr_spin
        # with spin, spin order: up up up dn dn dn
        elif ispin == 2:
            norbs = 2 * self.nwann
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:self.nwann, 0:self.nwann] = self.hr
            hr_spin[:, self.nwann:norbs, self.nwann:norbs] = self.hr
            return hr_spin
        # without spin
        else:
            return self.hr[:, :, :]


class KVec():
    """
    Define :math:`k` points in BZ, high symmetry line or uniform grid.

    Parameters
    ----------
    kpt_type: str
        The type of :math:`k` points, 'uni' or 'sym'.
    kbase: :math:`3 \\times 3` float array
        The basis vectors of the primitive reciprocal space.
    nkpt: int
        Number of :math:`k` points.
    kvec: float array
        The :math:`k` points.
    """

    def __init__(self, kpt_type='uni', kbase=None, nkpt=None, kvec=None):
        self.nkpt = nkpt
        self.kbase = np.array(kbase, dtype=np.float64)
        self.kvec = np.array(kvec, dtype=np.float64)
        self.kpt_type = kpt_type

    def set_base(self, kbase):
        """
        Set the basis of the primitive reciprocal.

        Parameters
        ----------
        kbase: :math:`3 \\times 3` float array
            The basis with respect to the global axis.
        """

        self.kbase = np.array(kbase, dtype=np.float64)

    def kvec_from_file(self, fname):
        """
        Read :math:`k` points from file.

        Parameters
        ----------
        fname: str
            File name.
        """

        tmp = []
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip().split()
                if line != []:
                    tmp.append(line)
            self.kvec = np.array(tmp, dtype=np.float64)
            self.nkpt = len(tmp)


class SymKVec(KVec):
    """
    Class for defining :math:`k` points in high symmetry line, derived from :class:`KVec`.

    Parameters
    ----------
    kbase: :math:`3 \\times 3` float array
        Basis of the primitive reciprocal lattice.
    hsymkpt: float array
        Starting and end :math:`k` points along high symmetry lines.
    klen: float array
        Length of segments of :math:`k` points line.
    """

    def __init__(self, kbase=None, hsymkpt=None, klen=None):
        self.klen = np.array(klen, dtype=np.float64)
        self.hsymkpt = np.array(hsymkpt, dtype=np.float64)
        KVec.__init__(self, 'sym', kbase)

    def get_klen(self):
        """
        Return length of :math:`k` points segments.
        """

        self.klen = np.zeros(self.nkpt, dtype=np.float64)
        self.klen[0] = 0.0
        prev_kpt = self.kvec[0]

        for i in range(1, self.nkpt):
            curr_kpt = self.kvec[i, :]
            tmp_kpt = curr_kpt - prev_kpt
            kx = np.dot(tmp_kpt, self.kbase[:, 0])
            ky = np.dot(tmp_kpt, self.kbase[:, 1])
            kz = np.dot(tmp_kpt, self.kbase[:, 2])
            self.klen[i] = self.klen[i - 1] + \
                np.sqrt(np.dot((kx, ky, kz), (kx, ky, kz)))
            prev_kpt = curr_kpt

    def from_hsymkpt(self, nkpt_per_path=20):
        """
        Given starting and end :math:`k` points of each segment,
        and the number of points per each segment,
        return the high symmetry :math:`k` points.

        Parameters
        ----------
        nkpt_per_path: int
            Number of :math:`k` points per each segment.
        """

        self.nkpt = nkpt_per_path * (len(self.hsymkpt) - 1)
        self.kvec = np.zeros((self.nkpt, 3), dtype=np.float64)
        for i in range(1, len(self.hsymkpt)):
            kpt_prev = self.hsymkpt[i - 1, :]
            kpt_curr = self.hsymkpt[i, :]
            for j in range(nkpt_per_path):
                ikpt = (i - 1) * nkpt_per_path + j
                self.kvec[ikpt, :] = (float(j) / float(nkpt_per_path - 1) *
                                      (kpt_curr - kpt_prev) + kpt_prev)

    def from_hsymkpt_uni(self, step):
        """
        Given a step, return high symmetry :math:`k` points.

        Parameters
        ----------
        step: float
            Step size.
        """

        kvec = []
        self.hsym_dis = np.zeros(len(self.hsymkpt), dtype=np.float64)
        self.hsym_dis[0] = 0.0
        for i in range(0, len(self.hsymkpt) - 1):
            kpt_prev = self.hsymkpt[i, :]
            kpt_curr = self.hsymkpt[i + 1, :]
            tmp = np.dot(self.kbase.transpose(), kpt_curr - kpt_prev)
            dis = np.sqrt(np.dot(tmp, tmp))
            self.hsym_dis[i + 1] = self.hsym_dis[i] + dis
            pts = np.arange(0, dis, step) / dis
            for ipt in pts:
                kvec.append(ipt * (kpt_curr - kpt_prev) + kpt_prev)
        self.kvec = np.array(kvec, dtype=np.float64)
        self.nkpt = len(self.kvec)


class UniKVec(KVec):
    """
    Class for defining uniform :math:`k` points grid, derived from :class:`KVec`.

    Parameters
    ----------
    grid: 3-elements tuple
        Three numbers defining a uniform grid, for example: :math:`11 \\times 11 \\times 11`.
    """

    def __init__(self, grid=None):
        self.grid = grid
        KVec.__init__(self, 'uni')

    def from_grid(self):
        """
        Return uniform :math:`k` points.
        """

        delta = 0.001
        nx, ny, nz = self.grid
        self.nkpt = nx * ny * nz
        self.kvec = np.zeros((self.nkpt, 3), dtype=np.float64)
        ikpt = 0
        for i in range(nx):
            if nx == 1:
                kx = 0.0
            else:
                kx = float(i) / float(nx)
            for j in range(ny):
                if ny == 1:
                    ky = 0.0
                else:
                    ky = float(j) / float(ny)
                for k in range(nz):
                    if nz == 1:
                        kz = 0.0
                    else:
                        kz = float(k) / float(nz)
                    ikpt = ikpt + 1
                    self.kvec[ikpt - 1, :] = kx + delta, ky + delta, kz + delta
