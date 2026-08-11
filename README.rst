===============================
edrixs
===============================

.. image:: https://img.shields.io/pypi/v/edrixs.svg
        :target: https://pypi.python.org/pypi/edrixs

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/NSLS-II/edrixs.git/master?urlpath=lab

EDRIXS is an open source toolkit for simulating XAS and RIXS spectra based on exact diagonalization of model Hamiltonians.
It was started as part of `COMSCOPE project <https://www.bnl.gov/comscope/software/comsuite.php/>`_ in the
Center for Computational Material Spectroscopy and Design, Brookhaven National Laboratory and is now maintained and
developed in collaboration between the `Condensed Matter Physics and Materials Science Division <https://www.bnl.gov/cmpmsd/>`_
and the `National Syncrotron Light Source II <https://www.bnl.gov/nsls2/>`_.

* Free software: GNU General Public License Version 3
* Documentation: https://edrixs.github.io/edrixs.
* Launch a `MyBinder Session <https://mybinder.org/v2/gh/edrixs/edrixs.git/master?urlpath=lab>`_ to try the code.

Features
--------

* ED solver
* XAS spectra
* RIXS spectra

How to cite
-----------
If you are using the EDRIXS code to do some studies and would like to publish your great works, it would be really appreciated if you can cite the following paper:

``EDRIXS: An open source toolkit for simulating spectra of resonant inelastic x-ray scattering, Y.L. Wang, G. Fabbris, M.P.M. Dean and G. Kotliar``, `Computer Physics Communications,243, 151 (2019) <https://doi.org/10.1016/j.cpc.2019.04.018>`_, `arXiv:1812.05735 <https://arxiv.org/abs/1812.05735/>`_.


Usage
-----
For Linux users we suggest installing with anaconda.

  .. code-block:: bash

     $ conda create --name edrixs_env python=3.10
     $ conda activate edrixs_env
     $ conda install -c conda-forge edrixs

 For Windows and macOS machines, we suggest using docker. See https://edrixs.github.io/edrixs for more details.
