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


This is the plan for some restructuring of the EDRIXS code

# Module organization
## solvers.py
Contains what the user will usually call to solve a particular problem. It should contain code that is independent of the method used to compute the spectra.
It is made up of functions
* `setup_<>` where `<>` is currently `1v1c`, `2v1c`, `siam`, etc. --- setup problem
* `get_H` --- build a many body Hamiltonian
* `ops` --- generate set of initial and intermediate state Hamiltonians and operators
* `ed`--- compute low-energy eigenenergies and eigenvectors. 
* `xas` --- generate final spectra
* `rixs` --- generate final spectra

`setup_<>` has negligible computational overhead and provides infrastructure to describe the problem, with many calls to the rest of the EDRIXS Python interface. `ops`, `ed`, `xas`, and `rixs` orchestrate calls into specific solvers. `get_H` will be called by `ops` and belongs in solvers.py since it is highly useful for the user. We plan to implement `get_H`, `ops`, `ed`, `xas`, `rixs` as thin wrapper functions. 
 
## <>_backend.py
Where <>= `scipy`, `petsc`, `numpy` (a dense matrix implementation) etc. A series of different modules re-implementing the underlying methods for different backends. _Each_ of these needs to re-implement the following functions.
* FockBasis_<> --- Class describing a FockBasis
* `get_H_<>` --- Build many-body Hamiltonian from basis and single-particle Hamiltonian. 
* `ed_<>`--- obtain low-energy eigenenergies and eigenvectors. 
* `xas_<>`
* `rixs_<>`

All of these will call into some underlying backend, which will depend on a few parameters. For uniformity, we should pass these as a dictionary `backend_kws`. Both `solvers.py` and `<>_backend.py` will implement Krylov methods, so there's no need to signal this in function names. 

The approach currently in `manybody_operator_csr.py` is backend-specific, so this algorithm should appear in scipy_backend.py as get_H_scipy(). This is because the Hamiltonian building process requires careful optimization that will probably be backend-dependent. 

## fock_basis.py
The `class FockBasis` currently in `manybody_operator_csr.py` belongs in fock_basis.py. I think manybody_operator_csr.py should be removed entirely and `fock_basis.py` should be reverted to its original form in the main branch. 

## solver_helpers.py
This should contain stuff that is back in independent and typically not called by the user. As mentioned below, I want to jettison the stuff handling parallelization, which is currently in the file and retain only functions like `_expand_broadening`
  

# Scheme for a typical calculation

## 1. Problem definition
Define the problem by calling `setup_<>` where `<>` is currently `1v1c`, `2v1c`, `siam`

```
emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat = setup_1v1c(...)
```

where `(...)` Defines the problem at hand _independent_ of which backend (scipy, petsc, etc.) is going to be called. 
`emat` (`umat`) is the two(four)-fermion Hamiltonian in the single-particle basis. Suffixes `_i`, `_n` denote the initial and intermediate states. `basis` describes (but does not build) the basis. `trans_mat` contains the three x,y,z transition operators from the core to valence states again in the single-particle basis. 

## 2. Build operators in the many-body basis. 

```
hmat_i, hmat_n, trans_ops = ops(emat_i, umat_i, basis_i, emat_n, umat_n, basis_n, trans_mat,
        *, backend=`scipy`)
```
`hmat_i` (`hmat_n`) is a many-body Hamiltonian for the initial (intermediate) states. trans_ops includes the x,y,z components of the many body absorption operator. 

Also within solvers, we should have a wrapper function

```
H = get_H(emat, umat, basis, backend=backend)
```

`ops` would essentially amount to three calls into `get_H`. Under the hood, `get_H` would call
for example `get_H_scipy` from `scipy_backend.py`.

From here, the matrices will be either scipy sparse or petsc. I would suggest that, from here, the functions call out to the appropriate backend based on type checking the input. Other options would be requiring passing a string specifying the solver or requiring the user to switch to different function names.



`get_H` would call out to get_H_scipy in 

## 3. Ground State Diagonalization

```
eval_all, evec_all = ed(hmat_i, ...)
```

Obtain low-energy eigenvalues and eigenvectors. 

## 4. XAS 
Not many notes needed; other than that, this should take a list of different polarizations. 

## 5. RIXS
This should take a 1d array of incident energies and polarizations. Note that if we want to parallelize this, we can easily pass in `[one single energy]` and `[one single polarization tuple]`

# Development plan
Work on the `petsc_edrixs` github branch. Raising PRs to `petsc_edrixs` as we have been doing. 

At least in the short term, I think we need to remove parallelization currently done with `from concurrent.futures import ProcessPoolExecutor`. There are a few reasons I say this, but perhaps the main one is that I want to defer the decision of how to do the parallelization until later in our process.

This can be done by making functions like
`rixs_scipy` run with simple for loops, but then moving threading into the test if needed. I.e. we want to expose a function without threading to the user. 

Some notes for tests could be written here...

Any particular problem you see for typing? 
     $ conda activate edrixs_env
     $ conda install -c conda-forge edrixs

 For Windows and macOS machines, we suggest using docker. See https://edrixs.github.io/edrixs for more details.
