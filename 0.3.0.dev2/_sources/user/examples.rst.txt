.. _examples:

Examples
========

There are two ways to learn EDRIXS from worked calculations. Start with the
pedagogical series to build up the underlying concepts and API one step at a
time, or go directly to the NiO examples for complete XAS and RIXS workflows.

Pedagogical series
------------------

The :doc:`pedagogical examples </auto_examples/index>` progress from exact
diagonalization, crystal fields, and Coulomb interactions to atomic and
Anderson impurity models. Each page includes the executable Python source and
its generated figures.

NiO XAS and RIXS
----------------

The following examples use the backend-independent model and solver APIs for
two descriptions of NiO:

.. toctree::
   :maxdepth: 1

   nio_l23_xas
   nio_siam_xas_rixs

The crystal-field calculation treats an ionic Ni :math:`3d^8` site. The
Anderson impurity calculation adds ligand orbitals and charge-transfer
configurations through :func:`~edrixs.models.model_siam`.
