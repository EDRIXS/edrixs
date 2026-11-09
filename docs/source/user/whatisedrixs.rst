.. _start-here:

Start here
==========

EDRIXS is an open source toolkit that uses exact diagonalization of model Hamiltonians to calculate X-ray
absorption (XAS) and resonant inelastic X-ray scattering (RIXS), as well as
energy levels, eigenstates, and observables of correlated-electron models.

It began as part of the `COMSCOPE project <https://www.bnl.gov/comscope/>`_
in the Center for Computational Material Spectroscopy and Design at Brookhaven
National Laboratory. It is maintained and developed in collaboration between
the `Condensed Matter Physics and Materials Science Division
<https://www.bnl.gov/cmpmsd/>`_ and the `National Synchrotron Light Source II
<https://www.bnl.gov/nsls2/>`_.

Choose a starting point for what you want to do:

* **Install EDRIXS:** :doc:`installation` covers released packages, source
  builds, and Colab, and links to the Docker and development-container guides.
  Use documentation for the version you install; the development documentation
  can describe interfaces that are absent from a released package.
* **Run a complete XAS or RIXS calculation:** start with
  :ref:`nio-l23-xas` for an ionic Ni model, then
  :ref:`nio-siam-xas-rixs` to include ligand hybridization and charge transfer.
  Both pages include executable Python code and the resulting spectra.
* **Understand the calculation workflow:** :ref:`introduction` explains how
  model parameters become Hamiltonians, eigenstates, and spectra.
* **Learn the physics or build a custom Hamiltonian:** the
  :doc:`/auto_examples/index` work through diagonalization, crystal fields,
  Coulomb interactions, and spectroscopy. They also illustrate calculations
  without core orbitals or photon transitions.
* **Check units, basis order, or spectral widths:** see :ref:`conventions`.
* **Choose a solver or improve convergence:** see :ref:`backend` and the
  backend settings linked there.
* **Look up a function:** see :ref:`reference`.

The scripts can run in a terminal or a notebook. You do not need to read every
pedagogical example before running the NiO calculations.
