************
Installation
************
For Linux users we suggest :ref:`installing with anaconda <AnacondaInstall>`. For Windows and macOS machines, we suggest using the :ref:`docker instructions <edrixsanddocker>`, which are relatively straightforward.  If desired, you can also compile the code from the source for Linux.


.. _AnacondaInstall:

Install and use edrixs via Anaconda
====================================================
A conda package has been built for Linux. To use edrixs via Anaconda, you need first to install `Anaconda <https://www.anaconda.com/distribution/>`_ in your system.
We recommend installing edrixs into a separate environment, for example, called ``edrixs_env``, together with any other packages you might want to use like this::

    conda create --name edrixs_env -c conda-forge python=3.10 edrixs matplotlib

We endeavor to keep the conda-forge release up to date, but note that these builds will usually not correspond to the latest version of edrixs, which is available in the `master branch of edrixs <https://github.com/NSLS-II/edrixs>`_.

edrixs will also run on `Google Colaboratory <https://research.google.com/colaboratory/>`_, but does not come installed as default. Installing it requires a you to install conda and then edrixs, which can be done by executing a cell::

    !pip install -q condacolab
    import condacolab
    condacolab.install()
    !conda install -c conda-forge edrixs

from within a notebook cell.

Requirements
============
Several tools and libraries are required to build and install edrixs,

   * Fortran compiler: gfortran and ifort are supported
   * MPI environment: openmpi and mpich are tested
   * MPI Fortran and C compilers: mpif90, mpicc
   * BLAS and LAPACK libraries: `OpenBLAS <https://github.com/xianyi/OpenBLAS/>`_ with gfortran and MKL with ifort
   * ARPACK library: `arpack-ng <https://github.com/opencollab/arpack-ng/>`_  with mpi enabled
   * Only Python3 is supported
   * numpy, scipy, sympy, matplotlib, sphinx, numpydoc
   * mpi4py with the same MPI implementation libraries (``openmpi`` or ``mpich``) as building edrixs

Build from source
=================
We will show how to build edrixs from source on Ubuntu Linux 20.04 and macOS Mojave (OSX 10.14) as examples.
We will use gcc, gfortran, openmpi and OpenBLAS in these examples.
Building edrixs on other versions of Linux or macOS, or with Intel's ifort+MKL will be similar.

Ubuntu Linux 20.04
------------------
Install compilers and tools::

    sudo apt-get update
    sudo apt-get install build-essential gfortran gcc
    sudo apt-get install git wget
    sudo apt-get install python3 libpython3-dev python3-pip python3-venv

Create and activate a python virtual environment for edrixs::

    python3 -m venv VIRTUAL_ENV
    source VIRTUAL_ENV/bin/activate

where ``VIRTUAL_ENV`` should be replaced by the directory where you wish to install edrixs.

Alternatively create and activate a conda environment for edrixs::

    conda create --name edrixs_env python=3.8
    conda activate edrixs_env

We will assume ``python`` and ``pip`` are pointing to the activated environment from now on.
Check we are using the expected python and pip::

    which python
    which pip
    python --version

Fetch the latest version of ``pip``::

    pip install --upgrade pip

openmpi, OpenBLAS, ARPACK can be installed by ``apt-get``, but their versions are old and may not work properly.
However, they can also be compiled from source easily. In the following, we will show both ways, but we always recommend to build newer ones from source.

openmpi can be installed by::

    sudo apt-get install libopenmpi-dev

or from newer version of source, for example v3.1.4::

    wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.4.tar.bz2
    tar -xjf openmpi-3.1.4.tar.bz2
    cd openmpi-3.1.4
    ./configure CC=gcc CXX=g++ FC=gfortran
    make
    sudo make install

the compiling process will take a while.

OpenBLAS can be installed by::

    sudo apt-get install libopenblas-dev

or from a newer version of source::

    wget https://github.com/xianyi/OpenBLAS/archive/v0.3.6.tar.gz
    tar -xzf v0.3.6.tar.gz
    cd OpenBLAS-0.3.6
    make CC=gcc FC=gfortran
    sudo make PREFIX=/usr/local install

ARPACK can be installed by::

    sudo apt-get install libarpack2-dev libparpack2-dev

or from a newer version of source::

    wget https://github.com/opencollab/arpack-ng/archive/3.6.3.tar.gz
    tar -xzf 3.6.3.tar.gz
    cd arpack-ng-3.6.3
    ./bootstrap
    ./configure --enable-mpi --with-blas="-L/usr/local/lib/ -lopenblas" FC=gfortran F77=gfortran MPIFC=mpif90 MPIF77=mpif90
    make
    sudo make install

mpi4py can be installed by::

    export MPICC=/usr/local/bin/mpicc
    sudo pip install --no-cache-dir mpi4py

or from source::

    wget https://github.com/mpi4py/mpi4py/archive/3.0.1.tar.gz
    tar xzf 3.0.1.tar.gz
    cd mpi4py-3.0.1

edit mpi.cfg to set MPI paths as following::

    [mpi]
    mpi_dir              = /usr/local
    mpicc                = %(mpi_dir)s/bin/mpicc
    mpicxx               = %(mpi_dir)s/bin/mpicxx
    include_dirs         = %(mpi_dir)s/include
    libraries            = mpi
    library_dirs         = %(mpi_dir)s/lib
    runtime_library_dirs = %(mpi_dir)s/lib

and comment all other contents. Then, build and install by::

    python setup.py build
    sudo pip install .

Check whether the MPI paths are correct by::

    python
    >>> import mpi4py
    >>> mpi4py.get_config()
    {'mpicc': '/usr/local/bin/mpicc',
     'mpicxx': '/usr/local/bin/mpicxx',
     'include_dirs': '/usr/local/include',
     'libraries': 'mpi',
     'library_dirs': '/usr/local/lib',
     'runtime_library_dirs': '/usr/local/lib'}

Now, we are ready to build edrixs::

    git clone https://github.com/NSLS-II/edrixs.git
    cd edrixs
    pip install -v .

Start to play with edrixs by::

    python
    >>> import edrixs
    >>> edrixs.some_functions(...)

or go to ``examples`` directory to run some examples::

    cd examples/more/ED/14orb
    ./get_inputs.py
    mpirun -np 2 ed.x
    mpirun -np 2 ./run_fedsolver.py
    cd ../../RIXS/LaNiO3_thin
    mpirun -np 2 ./run_rixs_fsolver.py

if no errors, the installation is successful.

macOS Mojave (OSX 10.14)
------------------------
Install newest Xcode through App store.

Use MacPorts
~~~~~~~~~~~~
Download and install `MacPorts <https://www.macports.org/install.php/>`_.
Update MacPorts by::

    sudo port -v selfupdate

Install gcc8, arpack, openblas and openmpi::

    sudo port -v install gcc8
    sudo port select gcc mp-gcc8
    sudo port -v install openmpi-default +gcc8
    sudo port -v install openblas +gcc8
    sudo port -v install arpack +openblas +openmpi
    sudo port select --set mpi openmpi-mp-fortran

Install Python, pip, numpy, scipy, sympy, matplotlib::

    sudo port -v install python37 py37-pip
    sudo port -v install py37-numpy +gcc8 +openblas
    sudo port -v install py37-scipy +gcc8 +openblas
    sudo port -v install py37-sympy
    sudo port -v install py37-matplotlib

**Notes:**

* DO NOT use pip to install numpy because it will use ``clang`` as default compiler, which has a strange bug when using ``f2py`` with ``mpif90`` compiler. If you cannot solve this issue by ``sudo port install py37-numpy +gcc8``, you can compile numpy from its source with ``gcc`` compiler. Always use gcc to compile numpy if you want to build it from source.

* You can also try ``gcc9`` if it is already available, but be sure to change all ``gcc8`` to ``gcc9`` in the above commands.

We will assume ``python`` pointing to ``python3.7`` and ``pip`` pointing to ``pip3.7`` from now on. If this is not the case, you can make links explicitly.
Check we are using the expected python and pip::

    which python
    python --version
    which pip
    pip --version

Add the following two lines into ``~/.bash_profile``::

    export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
    export PATH=/opt/local/Library/Frameworks/Python.framework/Versions/3.7/bin:$PATH

Close current terminal and open a new one.

Install mpi4py::

    export MPICC=/opt/local/bin/mpicc
    sudo pip install --no-cache-dir mpi4py

Please be sure to check whether the MPI paths of mpi4py are correct by::

    python
    >>> import mpi4py
    >>> mpi4py.get_config()
    {'mpicc': '/opt/local/bin/mpicc'}

Now, we are ready to build edrixs::

    git clone https://github.com/NSLS-II/edrixs.git
    cd edrixs
    make -C src F90=mpif90 LIBS="-L/opt/local/lib -lopenblas -lparpack -larpack"
    make -C src install
    python setup.py config_fc --f77exec=mpif90 --f90exec=mpif90 build_ext --libraries=openblas,parpack,arpack --library-dirs=/opt/local/lib
    sudo pip install .

You can add ``edrixs/bin`` to the environment variable ``PATH`` in ~/.bash_profile.

Go to ``examples`` directory to run some examples::

    cd examples/more/ED/14orb
    ./get_inputs.py
    mpirun -np 2 ../../../../src/ed.x
    mpirun -np 2 ./run_fedsolver.py
    cd ../../RIXS/LaNiO3_thin
    mpirun -np 2 ./run_rixs_fsolver.py

if no errors, the installation is successful.

All done, enjoy!

Use Homebrew
~~~~~~~~~~~~~
Install Homebrew::

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Add following line to ``~/.bash_profile``::

    export PATH="/usr/local/bin:$PATH"

Install gcc9::

    brew install gcc@9

Install openblas and arpack::

    brew install openblas
    brew install arpack

openmpi has been automatically installed when installing arpack.

Install python3.7::

    brew install python

We will assume ``python`` pointing to ``python3.7`` and ``pip`` pointing to ``pip3.7`` from now on. If this is not the case, you can make link explicitly.
Check we are using the expected python and pip::

    which python
    python --version
    which pip
    pip --version

Make links if gcc, g++ and gfortran are not pointing to gcc-9, g++-9, gfortran-9, for example::

    ln -s /usr/local/Cellar/gcc/9.1.0/bin/gcc-9 /usr/local/bin/gcc
    ln -s /usr/local/Cellar/gcc/9.1.0/bin/g++-9 /usr/local/bin/g++
    ln -s /usr/local/Cellar/gcc/9.1.0/bin/gfortran-9 /usr/local/bin/gfortran

DO NOT install numpy through ``pip`` because it uses ``clang`` as default compiler, which will cause problems.
We will build numpy from source with gcc::

    wget https://github.com/numpy/numpy/archive/v1.16.3.tar.gz
    tar xzf v1.16.3.tar.gz
    cd numpy-1.16.3
    export CC=gcc CXX=g++
    python setup.py build
    pip install .

You might need to do  ``brew install wget`` if it is not already installed.
If you have BLIS or MKL installed, you will need to tell numpy to compile with
openblas. Create a file in the numpy directory called site.cfg and put the
following text in it::

    [openblas]
    libraries = openblas
    library_dirs = /usr/local/Cellar/openblas/0.3.9/lib
    include_dirs = /usr/local/Cellar/openblas/0.3.9/include
    runtime_library_dirs = /usr/local/Cellar/openblas/0.3.9/lib

Now we are ready to install scipy, sympy, matplotlib::

    pip install scipy sympy matplotlib
    export MPICC=/usr/local/bin/mpicc
    pip install --no-cache-dir mpi4py

Please be sure to check whether the MPI paths of mpi4py are correct by::

    python
    >>> import mpi4py
    >>> mpi4py.get_config()
    {'mpicc': '/usr/local/bin/mpicc'}

Now, we are ready to build edrixs::

    git clone https://github.com/NSLS-II/edrixs.git
    cd edrixs
    make -C src F90=mpif90 LIBS="-L/usr/local/opt/openblas/lib -lopenblas -L/usr/local/lib -lparpack -larpack"
    make -C src install
    python setup.py config_fc --f77exec=mpif90 --f90exec=mpif90 build_ext --libraries=openblas,parpack,arpack --library-dirs=/usr/local/lib:/usr/local/opt/openblas/lib
    pip install .

You can add ``edrixs/bin`` to the environment variable ``PATH`` in ``~/.bash_profile``.

Go to ``examples`` directory to run some examples::

    cd examples/more/ED/14orb
    ./get_inputs.py
    mpirun -np 2 ../../../../src/ed.x
    mpirun -np 2 ./run_fedsolver.py
    cd ../../RIXS/LaNiO3_thin
    mpirun -np 2 ./run_rixs_fsolver.py

if no errors, the installation is successful.

All done, enjoy!

.. [#] To change your default pip you need to add a line to your ``~/.bashrc`` on linux or to your ``~/.bash_profile`` on macOS. This should be ``alias pip='/usr/bin/pip3'`` where the path is determined by calling ``which pip3`` from your terminal.
