.. _nio-siam-xas-rixs:

*******************************************
NiO Anderson impurity model XAS and RIXS
*******************************************

This example uses :func:`~edrixs.models.model_siam_2d1p` to calculate Ni
:math:`L_{2,3}`-edge XAS and RIXS for a NiO Anderson impurity model. The model
contains a correlated Ni :math:`3d` shell and one bath site representing ten
symmetry-adapted O :math:`2p` ligand spin-orbitals. Hybridization allows the
calculation to include ligand-to-metal charge-transfer configurations.
The distributed PETSc/SLEPc backend is used for the many-body operators and
solvers. The example runs on one process during the documentation build, or it
can be run in parallel with ``mpirun``.

The nominal impurity occupancy is ``nd=8``. The wrapper fixes one full
bath, so the initial valence sector contains 18 electrons; the intermediate
sector contains 19 valence electrons and five core electrons. The bath has
no Coulomb interactions.

The model parameters are specified before constructing the Hamiltonian:

* ``slater`` contains the initial three d-d integrals and the intermediate
  seven d-d and d-p integrals. Atomic multipolar integrals are scaled to 80%,
  and the monopole integrals are specified directly.
* ``Delta=4.7`` eV is the configuration-average energy cost of transferring
  one bath electron to the impurity, without crystal field, hybridization,
  or multipolar Coulomb interactions. It is not a bare shell-center splitting.
  The wrapper recovers the average interactions from ``slater`` and uses
  :func:`~edrixs.utils.CT_imp_bath` and
  :func:`~edrixs.utils.CT_imp_bath_core_hole` to set the shell centers.
* ``impurity_levels`` and ``bath_levels`` each contain five orbital energies
  ordered ``(dz2, dzx, dzy, dx2-y2, dxy)``. The mean of each set is subtracted
  before setting its center according to the charge transfer energy. Each value
  applies to both spins.
  Cubic splittings are 0.56 eV on the impurity and 1.44 eV on the bath.
* ``hyb`` gives five matching-orbital hoppings in that same order: 2.06 eV
  for the two eg orbitals and 1.21 eV for the three :math:`t_{2g}` orbitals.
* ``v_soc`` gives initial and intermediate impurity SOC; ``c_soc`` gives core
  SOC.
* ``ext_B`` applies an effective exchange field along [112], acting on spin.
  ``om_shift=857.6`` eV aligns the spectrum with the Ni L edge; the wrapper
  combines this with its internally derived core energy.

The first seven outputs from
:func:`~edrixs.models.model_siam_2d1p` pass to
:func:`~edrixs.solvers.get_ops`. Pass the eighth output, ``shift``, to
:func:`~edrixs.solvers.ed` to subtract the core shell-level energy during
diagonalization and restore it in the returned eigenvalues. The ED, XAS, and RIXS
steps select the SciPy backend. Temperature (300 K), photon geometry,
polarization, and lifetime broadening are supplied to the spectrum routines.

.. plot:: pyplots/nio_siam_xas_rixs.py
   :context: reset
   :include-source: True
