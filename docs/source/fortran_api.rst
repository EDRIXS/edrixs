.. _fortran_api:

###########
Fortran API
###########

Documentation for the EDRIXS Fortran back-end, extracted from source comments
via Doxygen and rendered here by Breathe.

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
