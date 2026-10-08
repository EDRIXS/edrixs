#!/usr/bin/env python
"""
Ground state analysis for NiO
================================================================================
This example follows the :ref:`sphx_glr_auto_examples_example_03_AIM_XAS.py`
example and considers the same model. This time we show how to analyze
the eigenvectors in terms of a
:math:`\\alpha |d^8L^{10}\\rangle + \\beta |d^9L^9\\rangle
+ \\gamma |d^{10}L^8\\rangle`
representation.
"""

################################################################################
# Hamiltonian
# ------------------------------------------------------------------------------
# We start by obtaining all eigenvectors of the previous Hamiltonian.

import numpy as np
import matplotlib.pyplot as plt
import edrixs
# sphinx_gallery_start_ignore
# Importing the script re-executes example 3. Silence its prints and discard any
# figures it leaves open so sphinx-gallery does not attach them to this page.
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    from example_03_AIM_XAS import emat_i, basis_i, hmat_i, norb_d, shift
plt.close('all')
# sphinx_gallery_end_ignore

backend = 'scipy'
eval_i, evec_i = edrixs.ed(hmat_i, shift=shift, num_evals=len(basis_i), backend=backend)
eval_i = eval_i - eval_i.min()

################################################################################
# Number of d electrons
# ------------------------------------------------------------------------------
# We can count the number of :math:`d` electrons by building a single-particle
# operator :code:`single_particle_Nd` with ones on the diagonal at the nickel
# :math:`d`-shell orbitals and zeros elsewhere. :func:`~edrixs.solvers.build_op`
# transforms this into the many-body operator :code:`O`, whose expectation value
# in each eigenvector is the :math:`d`-electron count of that state.

single_particle_Nd = np.zeros_like(emat_i)
single_particle_Nd[:norb_d, :norb_d] = np.eye(norb_d)

O = edrixs.build_op(single_particle_Nd, None, basis_i, backend=backend)

nd_expect = np.sum(evec_i.conj() * (O @ evec_i), axis=0).real

fig, ax = plt.subplots()

ax.set_xlabel('Energy (eV)')
ax.set_ylabel('$d$ electron count')
ax.plot(eval_i, nd_expect)

################################################################################
# Configuration weights alpha, beta, gamma
# ------------------------------------------------------------------------------
# Every Fock state has a definite :math:`d` occupation, so :math:`O` is diagonal
# in the Fock basis. :code:`O.diagonal()` therefore gives the :math:`d`-electron
# count of each Fock-basis component, in the same order as the amplitudes of the
# eigenvectors returned by :func:`~edrixs.solvers.ed`.
nd_per_component = np.rint(O.diagonal().real).astype(int)

amp2 = np.abs(evec_i)**2
alphas = amp2[nd_per_component == 8].sum(axis=0)
betas = amp2[nd_per_component == 9].sum(axis=0)
gammas = amp2[nd_per_component == 10].sum(axis=0)

print("Ground state\nalpha={:.3f}\tbeta={:.3f}\tgamma={:.3f}".format(
    alphas[0], betas[0], gammas[0]))

fig, ax = plt.subplots()
ax.plot(eval_i, alphas, label=r'$\alpha$ $d^8L^{10}$')
ax.plot(eval_i, betas, label=r'$\beta$ $d^9L^{9}$')
ax.plot(eval_i, gammas, label=r'$\gamma$ $d^{10}L^{8}$')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Population')
ax.legend()
