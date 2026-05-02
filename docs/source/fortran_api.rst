.. _fortran_api:

###########
Fortran API
###########

Documentation for the EDRIXS Fortran back-end, extracted from source comments
via Doxygen and rendered here by Breathe.

Overview of XAS and RIXS calculations
======================================

EDRIXS models spectroscopy on correlated electron systems using exact
diagonalisation (ED) of a many-body Hamiltonian written in a Fock-state basis.
Three distinct Hilbert spaces are involved, each with its own Hamiltonian and
Fock basis:

* **Initial space** (*i*): no core hole, described by hopping and Coulomb
  integrals read from ``hopping_i.in`` and ``coulomb_i.in``.
* **Intermediate space** (*n*): one core hole and one extra valence electron
  created by X-ray absorption (same total electron count as the initial space),
  described by ``hopping_n.in`` and ``coulomb_n.in``.
* **Final space** (*f*): used in RIXS only; no core hole, same total electron
  count as the initial space.  After the core hole is filled by an emitted
  photon the system is left in a (generally excited) valence configuration.
  The Hamiltonian has the same form as :math:`H_i`.

Step 1 — Exact diagonalisation (``ed.x`` / ``ed_driver``)
----------------------------------------------------------

The initial-state Hamiltonian :math:`H_i` is built in distributed CSR format
and diagonalised to find the lowest few eigenstates
:math:`|\psi_k\rangle` with energies :math:`E_k`.  Three solvers are
available (controlled by ``ed_solver``): full LAPACK ``ZHEEV`` for small
systems, a parallel Lanczos algorithm, or ARPACK.  The eigenvectors are
written to ``eigvec.k`` for use by the spectroscopy drivers.

Step 2 — Absorption transition operator
----------------------------------------

The dipole (or multipole) transition operator :math:`T` couples the initial
and intermediate spaces.  Its matrix elements are computed in
``build_transop_i`` by iterating over pairs of Fock states: for each
initial state the operator creates a core hole in one of the ``num_core_orbs``
core orbitals, searches the intermediate Fock basis for the resulting state
(with Jordan–Wigner sign), and stores the amplitude as a CSR entry.  The
result is a rectangular CSR matrix with ``ndim_n`` rows and ``ndim_i``
columns.

XAS (``xas.x`` / ``xas_driver``)
----------------------------------

The X-ray absorption spectrum is the imaginary part of the one-particle
Green's function

.. math::

   A(\omega) = -\frac{1}{\pi}\,\mathrm{Im}\,
   \langle\phi|({\omega - H_n + i\gamma})^{-1}|\phi\rangle,

where :math:`|\phi\rangle = T|\psi_k\rangle` is the intermediate state
produced by applying the transition operator to ground state :math:`k`, and
:math:`\gamma` is the core-hole lifetime broadening.

Rather than inverting :math:`\omega - H_n` explicitly, EDRIXS evaluates
:math:`A(\omega)` via the **Lanczos continued-fraction method**
(``build_krylov_mp``):

1. Run the Lanczos recurrence starting from the normalised vector
   :math:`|\phi\rangle / \||\phi\rangle\|` for ``nkryl`` steps, producing
   diagonal coefficients :math:`\alpha_j` and off-diagonal coefficients
   :math:`\beta_j`.
2. The Green's function is then expressed exactly as the continued fraction

   .. math::

      G(z) = \frac{\||\phi\rangle\|^2}
             {z - \alpha_1 - \dfrac{\beta_1^2}
             {z - \alpha_2 - \dfrac{\beta_2^2}{\ddots}}}

   which ``build_spectrum`` evaluates on a dense frequency mesh
   :math:`z = \omega + E_k + i\gamma` with no further matrix operations.

The Krylov data (:math:`\alpha`, :math:`\beta`, norm, :math:`E_k`) are saved
to ``xas_poles.k`` so that the spectrum can be regenerated at any resolution
without re-running the expensive Lanczos iteration.

RIXS (``rixs.x`` / ``rixs_driver``)
--------------------------------------

The RIXS cross-section is proportional to

.. math::

   I(\omega_\text{in}, \omega_\text{loss}) \propto
   \sum_f \bigl|\langle f|T_f
   (\omega_\text{in} - H_n + E_k + i\gamma_\text{in})^{-1}
   T_i|k\rangle\bigr|^2
   \,\delta(\omega_\text{loss} - E_f + E_k),

which requires the resolvent :math:`(H_n - \omega)^{-1}|\phi\rangle` at a
single complex energy :math:`\omega = \omega_\text{in} + E_k + i\gamma_\text{in}`.
This is obtained by solving the **complex-symmetric linear system**

.. math::

   (H_n - \omega)\,x = |\phi\rangle

with the parallel MINRES solver (``pminres_csr``).  The shifted Hamiltonian
stored in ``ham_csr`` already absorbs the sign convention used by MINRES.

Once :math:`x = (H_n - \omega)^{-1}|\phi\rangle` is known, the emission
transition operator :math:`T_f` (built by ``build_transop_f``) maps it into
the final-state space:

.. math::

   |\phi_f\rangle = T_f\,x.

The RIXS spectrum is then the XAS problem again: run ``build_krylov_mp``
on :math:`H_f` (same form as :math:`H_i`) starting from :math:`|\phi_f\rangle`,
write the poles to ``rixs_poles.k``, and evaluate the continued fraction
with ``build_spectrum`` for the desired energy-loss mesh.

MPI parallelism
----------------

All three spaces can be large (millions of basis states).  The Hamiltonian and
transition-operator matrices are distributed across MPI ranks in **row blocks**:
each rank owns a contiguous slice of rows and the diagonal column block.
Off-diagonal column blocks needed for the sparse matrix–vector product
(``pspmv_csr``) are communicated via non-blocking ``MPI_ISEND`` / ``MPI_RECV``
while the local diagonal multiply proceeds, overlapping communication and
computation.  Dot products and norms are reduced with ``MPI_ALLREDUCE``.

.. note::

   For the full Doxygen HTML output (including derived-type definitions from
   ``m_types`` and cross-referenced source browsing) run ``doxygen Doxyfile``
   from the repository root and open ``docs/doxygen/html/index.html``.

Modules
=======

m_constants
-----------

.. doxygennamespace:: m_constants
   :project: edrixs

m_control
---------

.. doxygennamespace:: m_control
   :project: edrixs

m_global
--------

.. doxygennamespace:: m_global
   :project: edrixs

m_lanczos
---------

.. doxygennamespace:: m_lanczos
   :project: edrixs

Source files
============

fock.f90
--------

.. doxygenfile:: fock.f90
   :project: edrixs

ham.f90
-------

.. doxygenfile:: ham.f90
   :project: edrixs

utils.f90
---------

.. doxygenfile:: utils.f90
   :project: edrixs

spmv.f90
--------

.. doxygenfile:: spmv.f90
   :project: edrixs

full_diag.f90
-------------

.. doxygenfile:: full_diag.f90
   :project: edrixs

linsys.f90
----------

.. doxygenfile:: linsys.f90
   :project: edrixs

io.f90
------

.. doxygenfile:: io.f90
   :project: edrixs

Solver drivers
==============

ed_driver.f90
-------------

.. doxygenfile:: ed_driver.f90
   :project: edrixs

xas_driver.f90
--------------

.. doxygenfile:: xas_driver.f90
   :project: edrixs

rixs_driver.f90
---------------

.. doxygenfile:: rixs_driver.f90
   :project: edrixs

opavg_driver.f90
----------------

.. doxygenfile:: opavg_driver.f90
   :project: edrixs
