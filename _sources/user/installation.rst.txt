************
Installation
************
For Linux and macOS, we recommend installation via Mamba. Windows users should use the :ref:`Docker instructions
<edrixsanddocker>`.


.. _AnacondaInstall:

Install the released package with Mamba on Linux and macOS
==========================================================

A prebuilt conda package is available for Linux and macOS. First install `Miniforge
<https://github.com/conda-forge/miniforge>`_, which provides Conda and Mamba,
then create a separate environment for edrixs::

    mamba create --name edrixs_env -c conda-forge --strict-channel-priority \
        edrixs matplotlib
    conda activate edrixs_env

edrixs will also run on `Google Colaboratory
<https://research.google.com/colaboratory/>`_, but is not installed there by
default. Install Conda and then edrixs from within a notebook cell::

    !pip install -q condacolab
    import condacolab
    condacolab.install()
    !conda install -c conda-forge edrixs

.. _MacOSInstall:
.. _SourceInstall:

Build the current source with Mamba on Linux or macOS
=====================================================

This is the recommended way to build the current source on Linux and macOS.
The same environment provides Python, OpenMPI, OpenBLAS, parallel ARPACK, and
gfortran, preventing incompatible system and conda libraries from being mixed.
On macOS, in particular, do not mix this environment with Homebrew or MacPorts
versions of these libraries, as that can cause linker or OpenMP runtime errors.

Install `Miniforge <https://github.com/conda-forge/miniforge>`_, then create
and activate the build environment::

    mamba create --name edrixs_env -c conda-forge --strict-channel-priority \
        python=3.14 "numpy>=2" scipy sympy matplotlib sphinx mpi4py \
        "arpack=*=mpi_openmpi*" openmpi gfortran \
        "libblas=*=*openblas" cmake ninja pip setuptools wheel
    conda activate edrixs_env

Clone and install edrixs from the repository root::

    git clone https://github.com/EDRIXS/edrixs.git
    cd edrixs
    python -m pip install --no-build-isolation --no-deps .

The ``--no-build-isolation`` option makes the extension use NumPy and the
compiler toolchain from the activated environment. The ``--no-deps`` option
prevents pip from replacing the compatible conda packages installed above.

Run these checks from outside the source directory so that Python imports the
installed package::

    cd ..
    python -c "from edrixs import fedrixs; print(fedrixs.__file__)"
    mpirun -np 2 python -c "from mpi4py import MPI; from edrixs import fedrixs; print(MPI.COMM_WORLD.rank)"

The first command should print the path to a file ending in ``.so``. The
second should print ranks ``0`` and ``1``. On macOS, the linker may print
harmless ``compact unwind`` warnings while compiling the gfortran objects.

Requirements
============
The Mamba command above installs compatible versions of all build and runtime
requirements. The supported versions and required components are:

   * Python 3.10 or newer; Python 3.14 is used in the recommended environment
   * NumPy 1.26 or newer at runtime; NumPy 2 or newer is used to build the
     extension for compatibility with both NumPy 1.26 and 2.x
   * SciPy, SymPy, Matplotlib, Sphinx, and numpydoc
   * CMake 3.17.3 or newer and Ninja
   * A Fortran compiler; the recommended environment currently uses gfortran
     16
   * An MPI environment with Fortran and C compilers; OpenMPI 5 and its
     ``mpif90`` and ``mpicc`` wrappers are used in the recommended environment
   * mpi4py 4 or newer, built with the same MPI implementation used to build
     edrixs
   * BLAS and LAPACK; OpenBLAS 0.3 is used in the recommended environment
   * `ARPACK-NG <https://github.com/opencollab/arpack-ng/>`_ 3.9 or newer,
     built with MPI support

Install with Docker on Windows
==============================

For Windows, we recommend using the maintained edrixs Docker image instead of
building the Fortran extension natively. See :ref:`edrixsanddocker` for the
image, Docker Compose configuration, and usage instructions.
