===============
Release History
===============

Unreleased
----------

New calculation interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Introduce a staged workflow that separates physical model construction,
  many-body operator construction, diagonalization, and spectral calculations:
  ``model_*`` → ``get_ops`` → ``ed`` → ``xas`` / ``rixs``.
* Add backend-independent ``model_1v1c``, ``model_2v1c``, ``model_siam``, and
  ``model_siam_2d1p`` constructors, including sparse Coulomb tensors and a
  charge-transfer parameterization for a d impurity, bath, and p core.
* Make SciPy the default backend, with sparse operator construction and Krylov
  methods for ED, XAS, and RIXS. Dense calculations and native Fortran/MPI
  solvers are also available through the staged interface.
* A PETSc-based backend is under development
* Add compact Fock-basis specifications with interchangeable combinadic and
  explicit basis construction. Operator construction uses Numba by default;
  ``use_numba=False`` selects Python construction.
* Validate ``backend_kws`` by backend and operation, with typo suggestions,
  documented defaults, and value constraints.

Numerical methods and utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Support an energy-reference shift for SciPy and dense diagonalization while
  returning eigenvalues in the original reference. Models return the shift
  alongside their operator inputs. The Fortran backend warns and continues
  when a nonzero shift is supplied.
* Reuse SciPy RIXS intermediate-state solves across outgoing polarization
  channels sharing the same incoming polarization.
* Add conversions between Racah parameters and d-shell Slater integrals.
* Expand tests for model construction, basis ordering, fermionic signs,
  backend consistency, and XAS/RIXS spectra.

Compatibility and migration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Deprecate the legacy model-specific Python and Fortran solver wrappers in
  favor of the staged interface. Existing wrappers emit deprecation warnings.
* Reject isotropic RIXS polarization; specify incoming and outgoing
  polarization channels explicitly. Isotropic XAS remains supported.


Documentation and packaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Update the examples for the staged workflow and document the SciPy, dense,
  and Fortran backend options and numerical-precision considerations.
* Clarify Slater-integral normalization and add versioned documentation and
  stricter documentation validation.
* Correct the CMake extension name used by packaging, automate Docker Hub
  release publishing, and update installation guidance and dependencies.

This page summarizes the user-visible changes in each tagged release.  See the
`GitHub releases page <https://github.com/EDRIXS/edrixs/releases>`_ for the
complete release notes and downloadable artifacts.

v0.1.3 (2026-09-08)
--------------------

* Updated the build system to support installation on Apple silicon.  This
  release does not change the physics functionality.
* Revised the charge-transfer examples and documentation and made minor
  example-maintenance improvements.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.1.3>`__ |
`Changes since v0.1.2 <https://github.com/EDRIXS/edrixs/compare/v0.1.2...v0.1.3>`__

v0.1.2 (2026-05-06)
--------------------

* Fixed several small issues in the solvers, examples, and documentation,
  including a missing spinor and outdated ``pkg_resources`` usage.
* Improved wheel generation, supported Python versions, development-container
  setup, installation guidance, and documentation publishing.
* Added and clarified examples for scattering geometry, NiO, and execution-time
  estimates.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.1.2>`__ |
`Changes since v0.1.1 <https://github.com/EDRIXS/edrixs/compare/v0.1.1...v0.1.2>`__

v0.1.1 (2024-06-05)
--------------------

* Fixed the ``end_indx`` handling in the Fortran expectation-value driver and
  corrected macOS compilation.
* Updated Docker packaging, installation documentation, and the publication
  list.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.1.1>`__ |
`Changes since v0.1.0 <https://github.com/EDRIXS/edrixs/compare/v0.1.0...v0.1.1>`__

v0.1.0 (2023-10-13)
--------------------

* Added a CMake-based build and improved macOS support.
* Added examples covering Hund's coupling and charge-transfer excitons, and
  refreshed the existing examples and Binder setup.
* Expanded the documentation and publication list.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.1.0>`__ |
`Changes since v0.0.8 <https://github.com/EDRIXS/edrixs/compare/v0.0.8...v0.1.0>`__

v0.0.8 (2022-10-28)
--------------------

* Fixed polarization normalization in the Python XAS solver.
* Added and revised RIXS, Hubbard-model, crystal-field, and transition examples,
  along with substantial documentation and continuous-integration updates.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.0.8>`__ |
`Changes since v0.0.7 <https://github.com/EDRIXS/edrixs/compare/v0.0.7...v0.0.8>`__

v0.0.7 (2020-09-07)
--------------------

* Improved the examples, documentation, and installation instructions.
* Added the missing Cu ``d10_2p5`` Slater parameters and fixed several small
  issues.

`Changes since v0.0.6 <https://github.com/EDRIXS/edrixs/compare/v0.0.6...v0.0.7>`__

v0.0.6 (2020-07-08)
--------------------

* Added three gallery examples, including a charge-transfer calculation.
* Added utilities that derive impurity and bath energy levels from the
  charge-transfer energy and ``U``, and fixed several bugs.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.0.6>`__ |
`Changes since v0.0.5 <https://github.com/EDRIXS/edrixs/compare/v0.0.5...v0.0.6>`__

v0.0.5 (2020-06-03)
--------------------

* Updated the documentation and examples, including expectation-value and
  Python-tips material.

`Release notes <https://github.com/EDRIXS/edrixs/releases/tag/v0.0.5>`__ |
`Changes since v0.0.4 <https://github.com/EDRIXS/edrixs/compare/v0.0.4...v0.0.5>`__

v0.0.4 (2019-07-01)
--------------------

* Corrected the source distribution and extension-linker configuration for
  packaging, including Conda packaging.

`Changes since v0.0.3 <https://github.com/EDRIXS/edrixs/compare/v0.0.3...v0.0.4>`__

v0.0.3 (2019-06-28)
--------------------

* Included ``requirements.txt`` in the source distribution.

`Changes since v0.0.2 <https://github.com/EDRIXS/edrixs/compare/v0.0.2...v0.0.3>`__

v0.0.2 (2019-06-27)
--------------------

* Added simpler APIs for angular-momentum and many-body operators.
* Corrected isotropic XAS calculations and improved Docker/Jupyter support.

`Changes since v0.0.1 <https://github.com/EDRIXS/edrixs/compare/v0.0.1...v0.0.2>`__

v0.0.1 (2019-05-15)
--------------------

* Initial tagged release of EDRIXS.

`Source at v0.0.1 <https://github.com/EDRIXS/edrixs/tree/v0.0.1>`__
