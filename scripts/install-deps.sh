#!/bin/bash

set -vxeuo pipefail

# The Python package is now pure Python. Native Fortran solvers are built and
# distributed separately, so this packaging/test environment does not need a
# Fortran compiler, MPI development libraries, CMake, or platform-specific
# wheel tagging.

# These packages are installed in the base environment but may be older
# versions. Explicitly upgrade them because they often create installation
# problems if out of date.
which python
python -VV

python -m pip install --upgrade pip setuptools wheel

# Generate source and universal Python wheel distributions.
python setup.py sdist bdist_wheel
ls -la dist/

# Install the wheel that was just built. The wheel is no longer tied to a
# CPython ABI tag such as cp311 because there is no compiled Python extension.
python -m pip install dist/edrixs-*.whl

# Install extra requirements for running tests and building docs.
python -m pip install -r requirements-dev.txt

# List the dependencies.
python -m pip list
