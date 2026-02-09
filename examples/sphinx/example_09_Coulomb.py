#!/usr/bin/env python
"""
Coulomb interactions
=====================================
In this example we provide more details on how Coulomb interactions are
implemented in multiplet calculations and EDRIXS in particular. We aim
to clarify the form of the matrices, how they are parametrized,
and how the breaking of spherical symmetry can switch on additional elements
that one might not anticipate. Our example is based on a :math:`d` atomic shell.
"""

################################################################################
# Create matrix
# ------------------------------------------------------------------------------
# The Coulomb interaction between two particles can be written as
#
#   .. math::
#     \begin{equation}
#     \hat{H} = \frac{1}{2}
#     \int d\mathbf{r} \int d\mathbf{r}^\prime
#     \Sigma_{\sigma, \sigma^\prime}
#     |\hat{\psi}^\sigma(\mathbf{r})|^2 \frac{e^2}{R}
#     |\hat{\psi}^{\sigma^\prime}(\mathbf{r^\prime})|^2,
#     \end{equation}
#
# where :math:`\hat{\psi}^\sigma(\mathbf{r})` is the electron wavefunction, with
# spin :math:`\sigma`, and :math:`R=|r-r^\prime|` is the electron separation.
# Solving our problem in this form is difficult due to the need to symmeterize
# the wavefunction to follow fermionic statistics.
# Using second quantization, we can use operators to impose the required
# particle exchange statistics and write the equation in terms of
# a tensor :math:`U`
#
#   .. math::
#     \begin{equation}
#     \hat{H} = \sum_{\alpha,\beta,\gamma,\delta,\sigma,\sigma^\prime}
#     U_{\alpha\sigma,\beta\sigma^\prime,\gamma\sigma^\prime,\delta\sigma}
#     \hat{f}^{\dagger}_{\alpha\sigma}
#     \hat{f}^{\dagger}_{\beta\sigma^\prime}
#     \hat{f}_{\gamma\sigma^\prime}\hat{f}_{\delta\sigma},
#     \end{equation}
#
# where :math:`\alpha`, :math:`\beta`, :math:`\gamma`, :math:`\delta` are
# orbital indices and :math:`\hat{f}^{\dagger}`
# (:math:`\hat{f}`) are the creation (anihilation) operators.
# For a :math:`d`-electron system, we have :math:`10` distinct spin-orbitals
# (:math:`5` orbitals each with :math:`2` spins), which makes matrix the
# :math:`10\times10\times10\times10` in total size.
# In EDRIXS the matrix can be created as follows:
import edrixs
import numpy as np
import scipy
import matplotlib.pyplot as plt
import itertools

F0, F2, F4 = 6.94, 14.7, 4.41
umat_chb = edrixs.get_umat_slater('d', F0, F2, F4)
################################################################################
# We stored this under variable :code:`umat_chb` where "cbh" stands for
# complex harmonic basis, which is the default basis in EDRIXS.

################################################################################
# Parameterizing interactions
# ------------------------------------------------------------------------------
# EDRIXS parameterizes the interactions in :math:`U` via Slater integral
# parameters :math:`F^{k}`. These relate to integrals of various spherical
# Harmonics as well as Clebsch-Gordon coefficients, Gaunt coefficients,
# and Wigner 3J symbols. Textbooks such as [1]_ can be used for further
# reference. If you are interested in the details of how
# EDRIXS does this (and you probably aren't) function :func:`.umat_slater`,
# constructs the required matrix via Gaunt coeficents from
# :func:`.get_gaunt`. Two alternative parameterizations are common.
# The first are the Racah parameters, which are
#
#   .. math::
#     \begin{eqnarray}
#     A &=& F^0 - \frac{49}{441} F^4 \\
#     B &=& \frac{1}{49}F^2 - \frac{5}{441}F^4 \\
#     C &=& \frac{35}{441}F^4.
#     \end{eqnarray}
#
# or an alternative form for the Slater integrals
#
#   .. math::
#     \begin{eqnarray}
#     F_0 &=& F^0  \\
#     F_2 &=& \frac{1}{49}F^2 \\
#     F_4 &=& \frac{1}{441}F^4,
#     \end{eqnarray}
#
# which involves different normalization parameters.

################################################################################
# Basis transform
# ------------------------------------------------------------------------------
# If we want to use the real harmonic basis, we can use a tensor
# transformation, which imposes the following orbital order
# :math:`3z^2-r^2, xz, yz, x^2-y^2, xy`, each of which involves
# :math:`\uparrow, \downarrow` spin pairs. Let's perform this transformation and
# store a list of these orbitals.
umat = edrixs.transform_utensor(umat_chb, edrixs.tmat_c2r('d', True))
orbitals = ['3z^2-r^2', 'xz', 'yz', 'x^2-y^2', 'xy']

################################################################################
# Interactions
# ------------------------------------------------------------------------------
# Tensor :math:`U` is a  series of matrix
# elements
#
#   .. math::
#     \begin{equation}
#     \langle\psi_{\gamma,\delta}^{\bar{\sigma},\bar{\sigma}^\prime}
#     |\hat{H}|
#     \psi_{\alpha,\beta}^{\sigma,\sigma^\prime}\rangle
#     \end{equation}
#
# the combination of which defines the energetic cost of pairwise
# electron-electron interactions between states :math:`\alpha,\sigma`
# and :math:`\beta,\sigma^\prime`. In EDRIXS we follow the convention of
# summing over all orbital pairs. Some other texts count each pair of
# indices only once. The matrix elements here will consequently
# be half the magnitude of those in other references.
# We can express the interactions in terms of
# the orbitals involved. It is common to distinguish "direct Coulomb" and
# "exchange" interactions. The former come from electrons in the same orbital
# and the later involve swapping orbital labels for electrons. We will use
# :math:`U_0` and :math:`J` as a shorthand for distinguishing these.
#
# Before we describe the different types of interactions, we note that since
# the Coulomb interaction is real, and due to the spin symmmetry properties
# of the process :math:`U` always obeys
#
#   .. math::
#     \begin{equation}
#     U_{\alpha\sigma,\beta\sigma^\prime,\gamma\sigma^\prime,\delta\sigma} =
#     U_{\beta\sigma,\alpha\sigma^\prime,\delta\sigma^\prime,\gamma\sigma} =
#     U_{\delta\sigma,\gamma\sigma^\prime,\beta\sigma^\prime,\alpha\sigma} =
#     U_{\gamma\sigma,\delta\sigma^\prime,\alpha\sigma^\prime,\beta\sigma}.
#     \end{equation}
#
#
# 1. Intra orbital
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The direct Coulomb energy cost to double-occupy an orbital comes from terms
# like :math:`U_{\alpha\sigma,\alpha\bar\sigma,\alpha\bar\sigma,\alpha\sigma}`.
# In this notation, we use :math:`\sigma^\prime` to denote that the matrix
# element is summed over all pairs and :math:`\bar{\sigma}` to denote sums
# over all opposite spin pairs. Due to rotational symmetry, all these
# elements are the same and equal to
#
#   .. math::
#     \begin{eqnarray}
#     U_0 &=& \frac{A}{2} + 2B + \frac{3C}{2}\\
#         &=& \frac{F_0}{2} + 2F_2 + 18F_4
#     \end{eqnarray}
#
# Let's print these to demonstrate where these live in the array
for i in range(0, 5):
    val = umat[i*2, i*2 + 1, i*2 + 1, i*2].real
    print(f"{orbitals[i]:<8} \t {val:.3f}")

################################################################################
# 2. Inter orbital Coulomb interactions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Direct Coulomb repulsion between different orbitals depends on terms like
# :math:`U_{\alpha\sigma,\beta\sigma^\prime,\beta\sigma^\prime,\alpha\sigma}`.
# Expresions for these parameters are provided in column :math:`U` in
# :ref:`table_2_orbital`. We can print the values from :code:`umat`
# like this:
for i, j in itertools.combinations(range(5), 2):
    val = umat[i*2, j*2 + 1, j*2 + 1, i*2].real
    print(f"{orbitals[i]:<8} \t {orbitals[j]:<8} \t {val:.3f}")

################################################################################
# 3. Inter-orbital exchange interactions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exchange terms exist with the form
# :math:`U_{\alpha\sigma,\beta\sigma^\prime,\alpha\sigma^\prime,\beta\sigma}`.
# Expresions for these parameters are provided in column :math:`J` of
# :ref:`table_2_orbital`. These come from terms like this in the matrix:
for i, j in itertools.combinations(range(5), 2):
    val = umat[i*2, j*2 + 1, i*2 + 1, j*2].real
    print(f"{orbitals[i]:<8} \t {orbitals[j]:<8} \t {val:.3f}")

################################################################################
# 4. Pair hopping term
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Terms that swap pairs of electrons exist as
# :math:`(1-\delta_{\sigma\sigma'})U_{\alpha\sigma,\alpha\bar\sigma,\beta\bar\sigma,\beta\sigma}`
# and depend on exchange interactions column :math:`J` from
# :ref:`table_2_orbital`
# and here in the matrix.
for i, j in itertools.combinations(range(5), 2):
    val = umat[i*2, i*2 + 1, j*2 + 1, j*2].real
    print(f"{orbitals[i]:<8} \t {orbitals[j]:<8} \t {val:.3f}")

################################################################################
# 5. Three orbital
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Another set of terms that one might not immediately anticipate involve three
# orbitals like
#
#   .. math::
#     \begin{equation}
#     U_{\alpha\sigma, \gamma\sigma', \beta\sigma', \gamma\sigma} \\
#     U_{\alpha\sigma, \gamma\sigma', \gamma\sigma', \beta\sigma} \\
#     (1-\delta_{\sigma\sigma'})
#     U_{\alpha\sigma, \beta\sigma', \gamma\sigma', \gamma\sigma}
#     \end{equation}
#
# for :math:`\alpha=3z^2-r^2, \beta=x^2-y^2, \gamma=xz/yz`.
# These are needed to maintain the rotational symmetry of the interations.
# See :ref:`table_3_orbital` for the expressions. We can print some of
# these via:
ijkl = [[0, 1, 3, 1],
        [0, 2, 3, 2],
        [1, 0, 3, 1],
        [1, 1, 3, 0],
        [2, 0, 3, 2],
        [2, 2, 3, 0]]

for i, j, k, l in ijkl:
    val = umat[i*2, j*2 + 1, k*2 + 1, l*2].real
    print(f"{orbitals[i]:<8} \t {orbitals[j]:<8} \t"
          f"{orbitals[k]:<8} \t {orbitals[l]:<8} \t {val:.3f}")

################################################################################
# 6. Four orbital
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Futher multi-orbital terms include
# :math:`U_{\alpha\sigma,\beta\sigma^\prime,\gamma\sigma^\prime,\delta\sigma}`.
# We can find these here in the matrix:
ijkl = [[0, 1, 2, 4],
        [0, 1, 4, 2],
        [0, 2, 1, 4],
        [0, 2, 4, 1],
        [0, 4, 1, 2],
        [0, 4, 2, 1],
        [3, 1, 4, 2],
        [3, 2, 4, 1],
        [3, 4, 1, 2],
        [3, 4, 2, 1]]

for i, j, k, l in ijkl:
    val = umat[i*2, j*2 + 1, k*2 + 1, l*2].real
    print(f"{orbitals[i]:<8} \t {orbitals[j]:<8} \t {orbitals[k]:<8}"
          f"\t {orbitals[l]:<8} \t {val:.3f}")

################################################################################
# Effects of multi-orbital terms
# ------------------------------------------------------------------------------
# To test the effects of the multi-orbital terms, let's plot the eigenenergy
# spectra with and without multi-orbital terms switched on for system with and
# without a cubic crystal field. We will use a :math:`d`-shell with two
# electrons.
ten_dqs = [0, 2, 4, 12]

def diagonalize(ten_dq, umat):
    emat = edrixs.cb_op(edrixs.cf_cubic_d(ten_dq),
                        edrixs.tmat_c2r('d', ispin=True))
    H = (edrixs.build_opers(4, umat, basis)
         + edrixs.build_opers(2, emat, basis))
    e, v = scipy.linalg.eigh(H)
    return e - e.min()

basis = edrixs.get_fock_bin_by_N(10, 2)
umat_no_multiorbital = np.copy(umat)
B = F2/49 - 5*F4/441
for val in [np.sqrt(3)*B/2, np.sqrt(3)*B, 3*B/2]:
    umat_no_multiorbital[(np.abs(umat)- val) < 1e-6] = 0

fig, axs = plt.subplots(1, len(ten_dqs), figsize=(8, 3))

for cind, (ax, ten_dq) in enumerate(zip(axs, ten_dqs)):
    ax.hlines(diagonalize(ten_dq, umat), xmin=0, xmax=1,
              label='on', color=f'C{cind}')
    ax.hlines(diagonalize(ten_dq, umat_no_multiorbital),
              xmin=1.5, xmax=2.5,
              label='off',
              linestyle=':', color=f'C{cind}')
    ax.set_title(f"$10D_q={ten_dq}$")
    ax.set_ylim([-.5, 20])
    ax.set_xticks([])
    ax.legend()

fig.suptitle("Eigenvalues with 3&4-orbital effects on/off")
fig.subplots_adjust(wspace=.3)
axs[0].set_ylabel('Eigenvalues (eV)')
fig.subplots_adjust(top=.8)
plt.show()

################################################################################
# On the left of the plot Coulomb interactions in spherical symmetry cause
# substantial mxing between :math:`t_{2g}` and :math:`e_{g}` orbitals in the
# eigenstates and 3 & 4 orbital orbital terms are crucial for obtaining the
# the right eigenenergies. As :math:`10D_q` get large, this mixing is switched
# off and the spectra start to become independent of whether the 3 & 4 orbital
# orbital terms are included or not.
#
#
#
# .. _table_2_orbital:
# .. table:: Table of 2 orbital interactions
#
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |Orbitals :math:`\alpha,\beta`|:math:`U_0` Racah | :math:`U_0` Slater    |:math:`J` Racah |:math:`J` Slater    |
#    +=============================+==================+=======================+================+====================+
#    |:math:`3z^2-r^2, xz`         |:math:`A/2+B+C/2` |:math:`F_0/2+F_2-12F_4`| :math:`B/2+C/2`|:math:`F_2/2+15F_4` |
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`3z^2-r^2, yz`         |:math:`A/2+B+C/2` |:math:`F_0/2+F_2-12F_4`| :math:`B/2+C/2`|:math:`F_2/2+15F_4` |
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`3z^2-r^2, x^2-y^2`    |:math:`A/2-2B+C/2`|:math:`F_0/2-2F_2+3F_4`|:math:`2B+C/2`  |:math:`2F_2+15F_4/2`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`3z^2-r^2, xy`         |:math:`A/2-2B+C/2`|:math:`F_0/2-2F_2+3F_4`|:math:`2B+C/2`  |:math:`2F_2+15F_4/2`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`xz, yz`               |:math:`A/2-B+C/2` |:math:`F_0/2-F_2-12F_4`|:math:`3B/2+C/2`|:math:`3F_2/2+10F_4`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`xz, x^2-y^2`          |:math:`A/2-B+C/2` |:math:`F_0/2-F_2-12F_4`|:math:`3B/2+C/2`|:math:`3F_2/2+10F_4`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`xz, xy`               |:math:`A/2-B+C/2` |:math:`F_0/2-F_2-12F_4`|:math:`3B/2+C/2`|:math:`3F_2/2+10F_4`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`yz, x^2-y^2`          |:math:`A/2-B+C/2` |:math:`F_0/2-F_2-12F_4`|:math:`3B/2+C/2`|:math:`3F_2/2+10F_4`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`yz, xy`               |:math:`A/2-B+C/2` |:math:`F_0/2-F_2-12F_4`|:math:`3B/2+C/2`|:math:`3F_2/2+10F_4`|
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#    |:math:`x^2-y^2, xy`          |:math:`A/2+2B+C/2`|:math:`F_0+4F_2-34F_4` | :math:`C/2`    |:math:`35F_4/2`     |
#    +-----------------------------+------------------+-----------------------+----------------+--------------------+
#
#
# .. _table_3_orbital:
# .. table:: Table of 3 orbital interactions
#
#    +-----------------------------+-------------+----------------------------------------------------+-----------------------------------------------------+
#    |Orbitals :math:`\alpha,\beta,\gamma,\delta`|:math:`\langle\alpha\beta|\gamma\delta\rangle` Racah|:math:`\langle\alpha\beta|\gamma\delta\rangle` Slater|
#    +=============================+=============+====================================================+=====================================================+
#    |:math:`3z^2-r^2, xz, x^2-y^2, xz`          | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, yz, x^2-y^2, yz`          | :math:`-\sqrt{3}B/2`                               | :math:`-\sqrt{3}F_2/2+5\sqrt{3}F_4/2`               |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`xz, 3z^2-r^2, x^2-y^2, xz`          | :math:`-\sqrt{3}B`                                 | :math:`-\sqrt{3}F_2+5\sqrt{3}F_4`                   |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`xz, xz, x^2-y^2, 3z^2-r^2`          | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`yz, 3z^2-r^2, x^2-y^2, yz`          | :math:`\sqrt{3}B`                                  | :math:`\sqrt{3}F_2-5\sqrt{3}F_4`                    |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`yz, yz, x^2-y^2, 3z^2-r^2`          | :math:`-\sqrt{3}B/2`                               | :math:`-\sqrt{3}F_2/2+5\sqrt{3}F_4/2`               |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#
#
# .. _table_4_orbital:
# .. table:: Table of 4 orbital interactions
#
#    +-----------------------------+-------------+----------------------------------------------------+-----------------------------------------------------+
#    |Orbitals :math:`\alpha,\beta,\gamma,\delta`|:math:`\langle\alpha\beta|\gamma\delta\rangle` Racah|:math:`\langle\alpha\beta|\gamma\delta\rangle` Slater|
#    +=============================+=============+====================================================+=====================================================+
#    |:math:`3z^2-r^2, xz, yz, xy`               | :math:`-\sqrt{3}B`                                 | :math:`-\sqrt{3}F_2+5\sqrt{3}F_4`                   |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, xz, xy, yz`               | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, yz, xz, xy`               | :math:`-\sqrt{3}B`                                 | :math:`-\sqrt{3}F_2+5\sqrt{3}F_4`                   |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, yz, xy, xz`               | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, xy, xz, yz`               | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`3z^2-r^2, xy, yz, xz`               | :math:`\sqrt{3}B/2`                                | :math:`\sqrt{3}F_2/2-5\sqrt{3}F_4/2`                |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`x^2-y^2 , xz, xy, yz`               | :math:`-3B/2`                                      | :math:`-3F_2/2+15F_4/2`                             |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`x^2-y^2 , yz, xy, xz`               | :math:`3B/2`                                       | :math:`3F_2/2-15F_4/2`                              |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`x^2-y^2 , xy, xz, yz`               | :math:`-3B/2`                                      | :math:`-3F_2/2+15F_4/2`                             |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#    |:math:`x^2-y^2 , xy, yz, xz`               | :math:`3B/2`                                       | :math:`3F_2/2-15F_4/2`                              |
#    +-------------------------------------------+----------------------------------------------------+-----------------------------------------------------+
#
#
# .. rubric:: Footnotes
#
# .. [1] MSugano S, Tanabe Y and Kamimura H. 1970. Multiplets of
#        Transition-Metal Ions in Crystals. Academic Press, New York and London.
