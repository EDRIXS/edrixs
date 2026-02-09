#!/usr/bin/env python
"""
Crystal fields
=====================================
This example explains how to implement crystal fields in edrixs.
"""
################################################################################
# We need to import these modules.
import edrixs
import numpy as np
import scipy

np.set_printoptions(precision=2, suppress=True, linewidth=90)

################################################################################
# Crystal field matrices
# ------------------------------------------------------------------------------
# Let us start by considering the common case of a :math:`d` atomic shell in a
# cubic crystal field. This is controlled by parameter :math:`10D_q` and is
# described in terms of a matrix which we will assign to :code:`cfmat`. edrixs
# can make this matrix via.
ten_dq = 10
cfmat = edrixs.angular_momentum.cf_cubic_d(ten_dq)

################################################################################
# Note this matrix is in a complex harmonic basis :math:`Y^m_l` where :math:`m`
# goes from :math:`-l,-l+1,...,l-1, l`. There is an up spin and a down
# spin for each :math:`Y^m_l`. This matrix is not diagonal in the complex
# harmonic basis, but it would be diagonal in the real harmonic basis
# :math:`d_{3z^2-r^2}, d_{xz}, d_{yz}, d_{x^2-y^2}, d_{xy}`.
# Let us diagonalize this matrix as a check and print out the energies
# and their degeneracies.

e, v = scipy.linalg.eigh(cfmat)
e = e.round(decimals=6)
unique_e = np.unique(e)
degeneracies = [sum(evalue == e) for evalue in unique_e]

print("E  \tDegeneracy")
for evalue, degenvalue in zip(unique_e, degeneracies):
    print("{:.1f}\t{:.0f}".format(evalue, degenvalue))
print("{} distinct energies".format(len(unique_e)))

################################################################################
# This makes sense! We see two different energies split by :math:`10D_q=10`. Let
# us look at the six columns corresponding to the lower energy eigenvalues.

print(v[:, :6].real)

################################################################################
# These are the set of so-called :math:`t_{2g}` orbitals, composed of
# :math:`Y^2_2, Y^{-2}_2, Y^{1}_2, Y^{-1}_2`. The rest of the eigenvectors
# (the last four) are
print(v[:, 6:].real)

################################################################################
# These are the set of so-called :math:`e_{g}` orbitals, composed of
# :math:`Y^2_2, Y^{-2}_2, Y^{0}_2`. We can use edrixs to prove that
# :code:`cfmat` would be diagonal in the real
# harmonic basis. An operator :math:`\hat{O}` can be transformed into an
# operator in another basis :math:`\hat{O}^{\prime}` using a unitary
# transformation matrix :math:`T` as
#
#    .. math::
#
#          \hat{O}^{\prime} = (T)^{\dagger} \hat{O} (T).
#
# This is computed as follows
cfmat_rhb = edrixs.cb_op(cfmat, edrixs.tmat_c2r('d', ispin=True))
print(cfmat_rhb.real)

################################################################################
# where :code:`edrixs.tmat_c2r('d', ispin=True)` is the transformation matrix.
# We needed to tell edrixs that we are working with a :math:`d`-shell and that it
# should include spin. We could also have transformed :code:`v` to see how these
# eignevectors are  composed of the real harmonic basis. We will see an example
# of this later.

################################################################################
# Crystal field on an atom
# ------------------------------------------------------------------------------
# To simulate the solid state, we need to combine the crystal field with Coulomb
# interactions. Let us choose an atomic model for Ni.
l = 2
norb = 10
noccu = 8
basis = edrixs.get_fock_bin_by_N(norb, noccu)
slater = edrixs.get_atom_data('Ni', '3d', noccu, edge='L3')['slater_i']

################################################################################
# Let us implement a tetragonal crystal field, for which we need to pass
# :code:`d1` the splitting of :math:`d_{yz}/d_{zx}` and :math:`d_{xy}` and
# :code:`d3` the splitting of :math:`d_{3z^2-r^2}` and :math:`d_{x^2-y^2}`.
ten_dq, d1, d3 = 2.5, 0.9, .2

################################################################################
# To determine the eigenvalues and eigenvectors we need to transform both our
# Coulomb matrix and our crystal field matrix into the same basis. See the
# example on exact diagonalization if needed. In this case, we put this
# procedure into a function, with the option to scale the Coulomb interactions.


def diagonlize(scaleU=1):
    umat = edrixs.get_umat_slater('d',
                                  slater[0][1]*scaleU,
                                  slater[1][1]*scaleU,
                                  slater[2][1]*scaleU)
    cfmat = edrixs.angular_momentum.cf_tetragonal_d(ten_dq=ten_dq, d1=d1, d3=d3)
    H = edrixs.build_opers(2, cfmat, basis) + edrixs.build_opers(4, umat, basis)
    e, v = scipy.linalg.eigh(H)
    e = e - np.min(e)  # define ground state as zero energy
    return e, v


################################################################################
# Let us look what happens when we run the function with the Coulomb
# interactions switched off and check the degeneracy of the output. Look at this
# python
# `string formatting tutorial <https://realpython.com/python-formatted-output/>`_
# if the code is confusing.
e, v = diagonlize(scaleU=0)
e = e.round(decimals=6)
unique_e = np.unique(e)
degeneracies = [sum(evalue == e) for evalue in unique_e]

print("E  \tDegeneracy")
for evalue, degenvalue in zip(unique_e, degeneracies):
    print("{:.1f}\t{:.0f}".format(evalue, degenvalue))
print("{} distinct energies".format(len(unique_e)))

################################################################################
# We see 10 distinct energies, which is the number of ways one can arrange
# two holes among 4 energy levels -- which makes sense as the tetragonal field
# involves four levels :math:`zx/zy, xy, 3z^2-r^2, x^2-y^2`. To see what is going
# on in more detail, we can also calculate the expectation
# values of the occupancy number of the orbitals
# :math:`3z^2-r^2, zx, zy, x^2-y^2, xy`.
# To create the operator, first write the matrix in the real harmonics basis
# :math:`|3z^2-r^2,\uparrow>`, :math:`|3z^2-r^2,\downarrow>`,
# :math:`|zx,\uparrow>`, :math:`|zx,\downarrow>`, etc.
# In this basis, they take a simple form: only the diagonal terms have element
# 1. We therefore make a 3D empty array and assign the diagonal as 1. Check
# out the
# `numpy indexing notes <https://numpy.org/doc/stable/reference/arrays.indexing.html>`_
# if needed.
nd_real_harmoic_basis = np.zeros((norb, norb, norb), dtype=complex)
indx = np.arange(norb)
nd_real_harmoic_basis[indx, indx, indx] = 1

################################################################################
# Recalling the necessity to put everything in the same basis, we transform
# into the complex harmonic basis and then transform into our Fock basis
nd_complex_harmoic_basis = edrixs.cb_op(nd_real_harmoic_basis,
                                        edrixs.tmat_r2c('d', True))
nd_op = edrixs.build_opers(2, nd_complex_harmoic_basis, basis)


################################################################################
# We apply the operator and print out as follows. Check the
# `numpy docs <https://numpy.org/doc/1.18/reference/generated/numpy.reshape.html>`_
# if the details of how the spin pairs have been added up is not immediately
# transparent.
nd_expt = np.array([edrixs.cb_op(nd_vec, v).diagonal().real for nd_vec in nd_op])

message = "{:>3s}" + "\t{:>6s}"*5
print(message.format(*"E 3z2-r2 zx zy x2-y2 xy".split(" ")))

message = "{:>3.1f}" + "\t{:>6.1f}"*5
for evalue, row in zip(e, nd_expt.T):
    spin_pairs = row.reshape(-1, 2).sum(1)
    print(message.format(evalue, *spin_pairs))

################################################################################
# The lowest energy state involves putting both holes in the :math:`x^2-y^2`
# orbital, which makes sense.  Now, let us redo the proceedure including Coulomb
# repulsion, which imposes an energy cost to putting multiple electrons in the
# same orbital.

e, v = diagonlize(scaleU=1)

nd_expt = np.array([edrixs.cb_op(nd_vec, v).diagonal().real for nd_vec in nd_op])

message = "{:>3s}" + "\t{:>6s}"*5
print(message.format(*"E 3z2-r2 zx zy x2-y2 xy".split(" ")))

message = "{:>3.1f}" + "\t{:>6.1f}"*5
for evalue, row in zip(e, nd_expt.T):
    spin_pairs = row.reshape(-1, 2).sum(1)
    print(message.format(evalue, *spin_pairs))

################################################################################
# Now the lowest energy state involves splitting the holes between the two
# highest energy :math:`x^2-y^2` and :math:`3z^2-r^2` orbitals. i.e. we have
# gone from low-spin to high spin. Working out the balance between these two
# states is often one of the first things one wants to determine upon the
# discovery of an interesting new material [1]_.

##############################################################################
#
# .. rubric:: Footnotes
#
# .. [1] M. Rossi et al., `arXiv:2011.00595  (2021) <https://arxiv.org/abs/2011.00595>`_.
