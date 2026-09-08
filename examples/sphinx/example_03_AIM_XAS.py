#!/usr/bin/env python
"""
Anderson impurity model for NiO XAS
================================================================================
Here we calculate the :math:`L`-edge XAS spectrum of an Anderson impurity model,
which is sometimes also called a charge-transfer multiplet model. This model
considers a set of correlated orbitals, often called the impurity or metal
states, that hybridize with a set of uncorrelated orbitals, often called the
ligands or bath states. Everyone's favorite test case for x-ray spectroscopic
calculations of the Anderson impurity model is NiO and we won't risk being
original! This means that our correlated states will
be the Ni :math:`3d` orbitals and the uncorrelated states come from O
:math:`2p` orbitals. The fact that we include these interactions means that our
spectrum can include processes where electrons transition from the bath to the
impurity, as such the Anderson Impurity Model is often more accurate than atomic
models, especially if the material has strong covalency.

When defining the bath states, it is useful to use the so-called
symmetry adapted linear combinations of orbitals as the basis. These states
take into account the symmetry relationships between the different bath atom
orbitals and the fact that there are bath orbital combinations that do not
interact with the impurity by symmetry. By doing this the problem can be
represented with fewer orbitals, which makes the calculation far more efficient.
The standard EDRIXS solver that we will use assumes that the bath states are
represented by an integer number of bath sites set by :code:`nbath`, each of
which hosts the same number of spin-orbitals as the impurity e.g. 10 for a
:math:`d`-electron material.

NiO has a rocksalt structure in which all Ni atoms are surrounded by six O
atoms. This NiO cluster used to simulate the
crystal would then contain 10 Ni :math:`3d` spin-orbitals and :math:`6`
spin-orbitals per O :math:`\\times 6` O atoms :math:`=36` oxygen
spin-orbitals. As explained by, for example, Maurits Haverkort
et al. in [1]_ symmetry allows us to represent the bath with 10 symmetry
adapted linear combinations of the different O :math:`p_x, p_y, p_z` states.
The crystal field and hopping parameters for
such a calculation can be obtained by post-processing DFT calculations. We will
use values for NiO from [1]_. If you use values from a paper the relevant
references should, of course, be cited.

"""
import edrixs
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# Number of electrons
# ------------------------------------------------------------------------------
# When formulating problems of this type, one usually thinks of a nominal
# valence for the impurity atom in this case :code:`nd = 8` and assume that the
# bath is full. The solver that we will
# use can simulate multiple bath sites. In our case we specify
# :code:`nbath  = 1` sites. Electrons will be able to transition from O to Ni
# during our calculation, but the total number of valence electrons
# :code:`v_noccu` will be conserved.
nd = 8
norb_d = 10
norb_bath = 10
nbath = 1
v_noccu  = nd + nbath*norb_d
shell_name = ('d', 'p') # valence and core shells for XAS calculation

################################################################################
# Coulomb interactions
# ------------------------------------------------------------------------------
# The atomic Coulomb interactions are usually initialized based on Hartree-Fock
# calculations from, for example,
# `Cowan's code <https://www.tcd.ie/Physics/people/Cormac.McGuinness/Cowan/>`_.
# edrixs has a database of these.
info  = edrixs.utils.get_atom_data('Ni', '3d', nd, edge='L3')

################################################################################
# The atomic values are typically scaled to account for screening in the solid.
# Here we use 80% scaling. Let's write these out in full, so that nothing is
# hidden. Values for :math:`U_{dd}` and :math:`U_{dp}` are those of Ref. [1]_
# obtained by comparing theory and experiment [2]_ [3]_.
scale_dd = 0.8
F2_dd = info['slater_i'][1][1] * scale_dd
F4_dd = info['slater_i'][2][1] * scale_dd
U_dd = 7.3
F0_dd = U_dd + edrixs.get_F0('d', F2_dd, F4_dd)

scale_dp = 0.8
F2_dp = info['slater_n'][4][1] * scale_dp
G1_dp = info['slater_n'][5][1] * scale_dp
G3_dp = info['slater_n'][6][1] * scale_dp
U_dp = 8.5
F0_dp = U_dp + edrixs.get_F0('dp', G1_dp, G3_dp)

slater = ([F0_dd, F2_dd, F4_dd],  # initial
          [F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp])  # with core hole

################################################################################
# Energy of the bath states
# ------------------------------------------------------------------------------
# In the notation used here, :math:`\Delta` sets the energy difference
# between the bath and impurity states. :math:`\Delta` is defined in the atomic
# limit without crystal field (i.e. in terms of the centers of the impurity and
# bath states before hybridization is considered) as the energy for a
# :math:`d^{n_d} \rightarrow d^{n_d + 1} \underline{L}` transition.
# Note that as electrons are moved one has to pay energy
# costs associated with the Coulomb interactions. The
# energy splitting between the bath and impurity is consequently not simply
# :math:`\Delta`. One must therefore determine the energies by solving
# a set of linear equations. See the :ref:`edrixs.utils functions <utils>`
# for details. We can call these functions to get the impurity energy
# :math:`E_d`, bath energy :math:`E_L`, impurity energy with a core hole
# :math:`E_{dc}`, bath energy with a core hole :math:`E_{Lc}` and the
# core hole energy :math:`E_p`. The initial ground state calculation is
# done in electron language leaving out the core shell.
Delta = 4.7
E_d, E_L = edrixs.CT_imp_bath(U_dd, Delta, nd)
################################################################################
# In the intermediate state, we include the core shell and the core hole
# potential. For this reason, the energies will shift differently to account
# for this.
E_dc, E_Lc, E_p = edrixs.CT_imp_bath_core_hole(U_dd, U_dp, Delta, nd)
message = ("E_d = {:.3f} eV\n"
           "E_L = {:.3f} eV\n"
           "E_dc = {:.3f} eV\n"
           "E_Lc = {:.3f} eV\n"
           "E_p = {:.3f} eV\n")
if __name__ == "__main__":
    print(message.format(E_d, E_L, E_dc, E_Lc, E_p))


################################################################################
# The spin-orbit coupling for the valence electrons in the ground state, the
# valence electrons with the core hole present, and for the core hole itself
# are initialized using the atomic values.
zeta_d_i = info['v_soc_i'][0]
zeta_d_n = info['v_soc_n'][0]
c_soc = info['c_soc']

################################################################################
# Build matrices describing interactions
# ------------------------------------------------------------------------------
# edrixs uses complex spherical harmonics as its default basis set. If we want to
# use another basis set, we need to pass a matrix to the solver, which transforms
# from complex spherical harmonics into the basis we use.
# The solver will use this matrix when implementing the Coulomb interactions
# using the :code:`slater` list of Coulomb parameters.
# Here it is easiest to
# use real harmonics. We make the complex harmonics to real harmonics transformation
# matrix via
trans_c2n = edrixs.tmat_c2r('d',True)

################################################################################
# The crystal field and SOC needs to be passed to the solver by constructing
# the impurity matrix in the real harmonic basis. For cubic symmetry, we need
# to set the energies of the orbitals along the
# diagonal of the matrix. These need to be in pairs as there are two
# spin-orbitals for each orbital energy. Python
# `list comprehension <https://realpython.com/list-comprehension-python/>`_
# and
# `numpy indexing <https://numpy.org/doc/stable/reference/arrays.indexing.html>`_
# are used here. See :ref:`sphx_glr_auto_examples_example_01_crystal_field.py`
# for more details if needed.
ten_dq = 0.56
CF = np.zeros((norb_d, norb_d), dtype=complex)
diagonal_indices = np.arange(norb_d)

orbital_energies = np.array([e for orbital_energy in
                             [+0.6 * ten_dq, # dz2
                              -0.4 * ten_dq, # dzx
                              -0.4 * ten_dq, # dzy
                              +0.6 * ten_dq, # dx2-y2
                              -0.4 * ten_dq] # dxy)
                             for e in [orbital_energy]*2])


CF[diagonal_indices, diagonal_indices] = orbital_energies

################################################################################
# The valence band SOC is constructed in the normal way and transformed into the
# real harmonic basis.
soc = edrixs.cb_op(edrixs.atom_hsoc('d', zeta_d_i), edrixs.tmat_c2r('d', True))

################################################################################
# The total impurity matrices for the ground and core-hole states are then
# the sum of crystal field and spin-orbit coupling. We further needed to apply
# an energy shift along the matrix diagonal, which we do using the
# :code:`np.eye` function which creates a diagonal matrix of ones.
E_d_mat = E_d*np.eye(norb_d)
E_dc_mat = E_dc*np.eye(norb_d)
imp_mat = CF + soc + E_d_mat
imp_mat_n = CF + soc + E_dc_mat

################################################################################
# The energy level of the bath(s) is described by a matrix where the row index
# denotes which bath and the column index denotes which orbital. Here we have
# only one bath, with 10 spin-orbitals. We initialize every entry to the bath
# energy :code:`E_L` and then split the levels according to :code:`ten_dq_bath`.
ten_dq_bath = 1.44
bath_level = np.full((nbath, norb_d), E_L, dtype=complex)
bath_level[0, :2] += ten_dq_bath*.6  # 3z2-r2
bath_level[0, 2:6] -= ten_dq_bath*.4  # zx/yz
bath_level[0, 6:8] += ten_dq_bath*.6  # x2-y2
bath_level[0, 8:] -= ten_dq_bath*.4  # xy
bath_level_n = np.full((nbath, norb_d), E_Lc, dtype=complex)
bath_level_n[0, :2] += ten_dq_bath*.6  # 3z2-r2
bath_level_n[0, 2:6] -= ten_dq_bath*.4  # zx/yz
bath_level_n[0, 6:8] += ten_dq_bath*.6  # x2-y2
bath_level_n[0, 8:] -= ten_dq_bath*.4  # xy

################################################################################
# The hybridization matrix describes the hopping between the bath
# and the impurity. This is called either :math:`V` or :math:`T` in the
# literature and matrix sign can either be positive or negative based.
# This is the same shape as the bath matrix. We take our
# values from Maurits Haverkort et al.'s DFT calculations [1]_.
Veg = 2.06
Vt2g = 1.21

hyb = np.zeros((nbath, norb_d), dtype=complex)
hyb[0, :2] = Veg  # 3z2-r2
hyb[0, 2:6] = Vt2g  # zx/yz
hyb[0, 6:8] = Veg  # x2-y2
hyb[0, 8:] = Vt2g  # xy

################################################################################
# We now need to define the parameters describing the XAS. X-ray polarization
# can be linear, circular or isotropic (appropriate for a powder).
poltype_xas = [('isotropic', 0)]
################################################################################
# edrixs uses the temperature in Kelvin to work out the population of the low-lying
# states via a Boltzmann distribution.
temperature = 300
################################################################################
# The x-ray beam is specified by the incident angle and azimuthal angle in radians
thin = 0 / 180.0 * np.pi
phi = 0.0
################################################################################
# these are with respect to the crystal field :math:`z` and :math:`x` axes
# written above. (That is, unless you specify the :code:`loc_axis` parameter
# described in the :code:`edrixs.model_siam` function documentation.)

################################################################################
# The spectrum in the raw calculation is offset by the energy involved with the
# core hole state, which is roughly :math:`5 E_p`, so we offset the spectrum by
# this and use :code:`om_shift` as an adjustable parameters for comparing
# theory to experiment. We also use this to specify :code:`ominc_xas`
# the range we want to compute the spectrum over. The core hole lifetime
# broadening also needs to be set via :code:`gamma_c`.
om_shift = 857.6
c_level = -om_shift - 5*E_p
ominc_xas = om_shift + np.linspace(-15, 25, 1000)

################################################################################
# The final state broadening is specified in terms of half-width at half-maximum
# You can either pass a constant value or an array the same size as
# :code:`om_shift` with varying values to simulate, for example, different state
# lifetimes for higher energy states.
gamma_c = np.full(ominc_xas.shape, 0.48/2)

################################################################################
# Magnetic field is a three-component vector in eV specified with respect to the
# same local axis as the x-ray beam. :code:`on_which = 'both'` applies the
# operator to the total spin plus orbital angular momentum, as is appropriate
# for a physical external magnetic field. Passing :code:`on_which = 'spin'`
# instead applies it to spin only, which is a convenient way to impose a
# magnetic order direction on the sample. The Bohr magneton
# :math:`\mu_B = 5.7883818012\times 10^{-5}` eV/T is useful for converting
# a physical field strength. Here we mimic magnetic order by applying a small
# spin field along :math:`z`.
ext_B = np.array([0.00, 0.00, 0.12])
on_which = 'spin'

################################################################################
# Build the model
# ------------------------------------------------------------------------------
# The staged EDRIXS interface separates the physical model from its numerical
# representation. :func:`~edrixs.model_siam` collects the one-body matrices,
# Coulomb parameters and Fock-basis specifications for the impurity plus bath
# problem. With :code:`siam_type=0` the model is assembled from :code:`imp_mat`,
# :code:`bath_level` and :code:`hyb` (and their core-hole counterparts). The
# external magnetic field is applied via :code:`ext_B` acting
# :code:`on_which='spin'`. The returned :code:`i` quantities describe the
# initial and final states without a core hole, while the :code:`n` quantities
# describe the intermediate state with a core hole. :code:`trans_mat` holds the
# Cartesian dipole matrices.
out = edrixs.model_siam(
    shell_name, nbath, siam_type=0, v_noccu=v_noccu,
    c_level=c_level, c_soc=c_soc, trans_c2n=trans_c2n,
    imp_mat=imp_mat, imp_mat_n=imp_mat_n,
    bath_level=bath_level, bath_level_n=bath_level_n, hyb=hyb,
    slater=slater, ext_B=ext_B, on_which=on_which,
)
emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat = out

################################################################################
# Diagonalization
# ------------------------------------------------------------------------------
# :func:`~edrixs.get_ops` converts these backend-independent ingredients into
# many-body initial/final and intermediate Hamiltonians, plus the dipole
# operators that map the initial Fock space to the intermediate one. Here we
# select the PETSc backend, which stores the Hamiltonians as distributed
# sparse matrices and uses SLEPc to extract the lowest eigenpairs. This scales
# to the larger Hilbert spaces typical of Anderson impurity models and runs in
# parallel if the script is launched with::
#
#        mpirun -n <number of processors> python example_03_AIM_XAS.py
#
# Running it as a plain :code:`python` script also works, it is just slower.
backend = 'petsc'
hmat_i, hmat_n, trans_ops = edrixs.get_ops(
    emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat,
    backend=backend,
)

################################################################################
# :func:`~edrixs.ed` obtains the retained low-energy eigenpairs of the
# Hamiltonian without a core hole. Here :code:`num_evals=3` states are
# thermally populated at the temperature of interest.
eval_i, evec_i = edrixs.ed(hmat_i, num_evals=3, backend=backend)

################################################################################
# Compute XAS
# ------------------------------------------------------------------------------
# The spectrum is built by applying the dipole operators to the thermally
# populated initial states and propagating in the intermediate-state
# Hamiltonian. EDRIXS weights the retained eigenstates by their Boltzmann
# factors at :code:`temperature`. Because our model includes hybridization,
# the spectrum captures charge-transfer processes in which electrons move
# between the O bath and the Ni impurity.
xas = edrixs.xas(
    eval_i, evec_i, hmat_n, trans_ops, ominc_xas,
    gamma_c=gamma_c, thin=thin, phi=phi, pol_type=poltype_xas,
    temperature=temperature, backend=backend,
)

################################################################################
# Let's plot the data and save it just in case. The returned array has shape
# :code:`(len(ominc_xas), len(poltype_xas))`.
if __name__ == "__main__":
    fig, ax = plt.subplots()

    ax.plot(ominc_xas, xas[:, 0])
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('XAS intensity')
    ax.set_title('Anderson impurity model for NiO')
    plt.show()

np.savetxt('xas.dat', np.column_stack((ominc_xas, xas)))

##############################################################################
#
# .. rubric:: Footnotes
#
# .. [1] Maurits Haverkort et al
#        `Phys. Rev. B 85, 165113 (2012) <https://doi.org/10.1103/PhysRevB.85.165113>`_.
# .. [2] A. E. Bocquet et al.,
#        `Phys. Rev. B 53, 1161 (1996) <https://doi.org/10.1103/PhysRevB.53.1161>`_.
# .. [3] Arata Tanaka, and Takeo Jo,
#        `J. Phys. Soc. Jpn. 63, 2788-2807(1994) <https://doi.org/10.1143/JPSJ.63.2788>`_.
