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

This is the recommended way to build the current source on
Linux and macOS. See the :ref:`macOS compiler note <macos-sdk-note>` if the
standard build fails on a recent macOS release.

Install `Miniforge <https://github.com/conda-forge/miniforge>`_, then create
and activate the build environment::

    mamba create --name edrixs_env -c conda-forge --strict-channel-priority \
        python=3.14 "numpy>=2" scipy sympy matplotlib mpi4py \
        "arpack=*=mpi_openmpi*" openmpi gfortran \
        "libblas=*=*openblas" cmake ninja pip setuptools wheel \
        sphinx ipython numpydoc pillow sphinx-copybutton sphinx-gallery \
        sphinx_rtd_theme pytest
    conda activate edrixs_env

Clone and install edrixs from the repository root::

    git clone https://github.com/EDRIXS/edrixs.git
    cd edrixs
    python -m pip install --no-build-isolation --no-deps .

The ``--no-build-isolation`` option makes the extension use NumPy and the
compiler toolchain from the activated environment. The ``--no-deps`` option
prevents pip from replacing the compatible conda packages installed above.

Requirements
============
The Mamba command above installs compatible versions of all build and runtime
requirements. The supported versions and required components are:

   * Python 3.10 or newer; Python 3.14 is used in the recommended environment
   * NumPy 1.26 or newer at runtime; NumPy 2 or newer is used to build the
     extension for compatibility with both NumPy 1.26 and 2.x
   * SciPy, SymPy, and Matplotlib
   * Sphinx, IPython, numpydoc, Pillow, sphinx-copybutton, sphinx-gallery, and
     sphinx-rtd-theme for building the documentation
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

.. _macos-sdk-note:

.. rubric:: macOS compiler note

If compilation reports that ``libSystem.tbd`` is malformed or contains an
unknown architecture such as ``arm64e.x1``, use Apple Clang and Apple's
SDK-compatible linker while retaining gfortran, OpenMPI, ARPACK, and BLAS from
the active Conda environment::

    xcode-select -p || xcode-select --install
    export PATH="$CONDA_PREFIX/bin:/usr/bin:/bin"
    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export FC="$CONDA_PREFIX/bin/gfortran"
    export F77="$FC"
    export OMPI_CC="$CC"
    export OMPI_CXX="$CXX"
    export OMPI_FC="$FC"
    apple_bin="$(dirname "$(xcrun --find ld)")"
    export FFLAGS="-B${apple_bin}"
    export FCFLAGS="$FFLAGS"
    export LDFLAGS="-B${apple_bin}"
    python -m pip install \
        --no-build-isolation --no-deps .

When retrying after changing compilers or environments, first remove the
generated CMake cache with ``cmake -E remove_directory build``.
