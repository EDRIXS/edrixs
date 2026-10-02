.. _conventions:

***********
Conventions
***********

This page summarizes the conventions used in EDRIXS: units and naming,
orbital orderings and single-particle bases, Fock-state encoding, and the
geometry, polarization, and broadening used to calculate spectra. The
:ref:`introduction` describes the interface, and the
:ref:`pedagogical examples <examples>` explore the underlying physics in
more detail.

Energies and naming
===================

* **Electron and hole language.** EDRIXS uses electron language by default.
  Using hole language involves inverting the energy levels in ``emat`` and
  taking the conjugate transpose of the absorption and emission operators.
* **Units.** All energies are in **eV**.
* **Initial and intermediate states.** Throughout EDRIXS, an ``_i`` suffix
  labels the initial and final states, which have no core hole, while ``_n``
  labels the intermediate state, which has one. Thus, ``emat_i``, ``umat_i``,
  ``basis_i``, and ``hmat_i`` describe the problem without a core hole, while
  ``emat_n``, ``umat_n``, ``basis_n``, and ``hmat_n`` describe the problem
  with a core hole.
* **Coulomb interactions.** These interactions are parameterized using
  Slater integrals. EDRIXS provides Hartree-Fock values for these integrals
  through :func:`~edrixs.utils.get_atom_data` and conversions between Slater
  integrals and the common :math:`U`, :math:`J`, :math:`U_d`, and :math:`J_H`
  parameterizations through functions such as
  :func:`~edrixs.utils.UdJH_to_F0F2F4` and
  :func:`~edrixs.utils.F0F2F4_to_UdJH`. The Hartree-Fock values are usually
  **scaled down** to 70-90% of their atomic values to approximate screening in
  the solid.
* **Core-level energies.** The absolute energy of a core level is not defined
  by the calculation. The resonance position is set through
  ``shell_level`` (or ``c_level``) together with an offset chosen to match
  experiment.

Core orbitals and energy shifts
===============================

EDRIXS model constructors include the core orbitals in the initial Fock space.
Because these orbitals are fully occupied, including them changes the state
labels without increasing the number of states. This convention avoids the
energy corrections needed when filled core orbitals are omitted and simplifies
the construction of transition operators.

Model constructors also return an energy offset, ``shift``. Passing it to
:func:`~edrixs.solvers.ed` with the dense, SciPy, or PETSc backend subtracts
the offset from the Hamiltonian diagonal before diagonalization and restores it to the
returned eigenvalues. This preserves the physics while reducing numerical
errors caused by subtracting large, nearly equal numbers.

Default orbital ordering and single-particle bases
==================================================

Unless specified otherwise, EDRIXS uses these orderings.

* **Spin:** interleaved, ``up, dn, up, dn, ...``.
* **Complex spherical harmonics** :math:`Y_l^m`: ``m = -l, -l+1, ..., l-1, l``.
* **Real spherical harmonics** (the Wannier90 ordering):

  - ``p``:  :math:`p_x, p_y, p_z`
  - ``d``:  :math:`d_{3z^2-r^2}, d_{xz}, d_{yz}, d_{x^2-y^2}, d_{xy}`
  - ``t2g``:  :math:`d_{xz}, d_{yz}, d_{xy}`
  - ``f``:  :math:`f_{z^3}, f_{xz^2}, f_{yz^2}, f_{z(x^2-y^2)}, f_{xyz},
    f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}`

* :math:`\lvert j, j_z \rangle` **basis** (SOC diagonal): the
  :math:`j = l-1/2` block first, then the :math:`j = l+1/2` block, each ordered
  :math:`-j, -j+1, ..., j`.

The **default single-particle basis used to define the Fock basis** is:

* complex spherical harmonics for ``p``, ``d``, ``t2g``, and ``f`` (``p`` and
  ``t2g`` share the same complex-harmonic basis);
* the :math:`\lvert j, j_z \rangle` basis for ``p12``, ``p32``, ``d32``,
  ``d52``, ``f52``, and ``f72``.

Helper functions that return matrices or tensors -- ``get_umat_slater``,
``get_trans_oper``, ``cf_cubic_d``, ``cf_tetragonal_d``, and
``cf_trigonal_t2g``, for example -- return them **in this default basis**.

.. important::

   You may express :code:`emat` and :code:`umat` in any single-particle basis,
   but everything entering the Hamiltonian must be expressed in that
   *same* basis. Transform one-body matrices with
   :func:`~edrixs.basis_transform.cb_op` and Coulomb tensors with
   :func:`~edrixs.basis_transform.transform_utensor`. The recommended
   practice is to keep the default basis and only transform the extra matrices
   you supply (for example, a crystal-field matrix) into it. ``cb_op`` applies

   .. math::

      \hat{O}^{\prime} = T^{\dagger}\, \hat{O}\, T ,

   Use ``T = edrixs.tmat_c2r('d', ispin=True)`` to transform from complex to
   real harmonics, or ``edrixs.tmat_r2c`` for the reverse transformation.

Fock basis
==========

* **Valence orbitals come first.** When a problem has both valence and core
  electrons, all valence spin-orbitals are indexed before all core
  spin-orbitals.
* **Encoding.** A Fock state is a string of 1s (occupied) and 0s (empty) over
  the single-particle spin-orbitals, stored as an integer for efficiency.
  Orbital zero occupies the least significant bit:
  orbital ``i`` has mask ``1 << i``. Occupation lists remain in increasing
  orbital order; ordinary binary integer displays place orbital zero at the
  right. Basis row ordering is unchanged.
* **Specification.** A basis is described compactly by
  :class:`~edrixs.fock_basis.FockBasisSpec`, built from ``(norb, nocc)`` pairs
  -- one per shell, each with fixed occupancy. For a single shell,
  :func:`~edrixs.fock_basis.get_fock_basis_int` returns a realized
  integer-encoded basis directly. :func:`~edrixs.solvers.build_op` and
  :func:`~edrixs.solvers.get_ops` accept either the specification or a realized
  basis.

Spectral conventions
====================

* **Geometry** follows Figure 1 of Y. Wang *et al.*,
  `Comput. Phys. Commun. 243, 151 (2019)
  <https://doi.org/10.1016/j.cpc.2019.04.018>`_. The incident and scattered
  beam angles are set by ``thin`` and ``thout``, and the azimuth by ``phi``.
  All angles are in radians and relative to the sample axes. By default, the
  crystal-field :math:`x, y, z` axes coincide with the lab frame; use
  ``loc_axis`` (model) or ``scatter_axis`` (solver) to change this.
* **Polarization** is given as a list of channels. For XAS, each entry is
  ``(kind, angle)`` with ``kind`` one of ``'linear'``, ``'left'``, ``'right'``
  or ``'isotropic'`` (the last for powders).  For RIXS each entry is a 4-tuple
  ``(in_kind, in_angle, out_kind, out_angle)``, where each kind is ``'linear'``,
  ``'left'``, or ``'right'``; sum the two outgoing channels when the experiment
  does not resolve emitted polarization.  ``'isotropic'`` is not supported for
  RIXS and is rejected with a warning.
* **Broadening** is a Lorentzian half width at half maximum.  ``gamma_c`` is the
  inverse core-hole lifetime (it dominates XAS and the incident-energy axis of
  RIXS); ``gamma_f`` is the final-state or resolution width on the RIXS
  energy-loss axis. ``gamma_c`` may be a scalar or an array with the same
  shape as ``ominc``; ``gamma_f`` may be a scalar or an array with the same
  shape as ``eloss``.
* **Temperature** (in K) sets Boltzmann weights over the retained low-energy
  initial states, so ``num_evals`` in ``ed`` must be large enough to cover all
  thermally populated states.
* Returned array shapes: ``xas`` has shape ``(len(ominc), len(pol_type))``;
  ``rixs`` has shape ``(len(ominc), len(eloss), len(pol_type))``.
