#!/bin/bash

set -vxeuo pipefail

platform="$(uname)"

if [ "${platform}" == "Linux" ]; then
    # sudo apt-get update -y
    # sudo apt-get install -y \
    #     gfortran \
    #     openmpi-bin \
    #     openmpi-doc \
    #     libopenmpi-dev \
    #     libopenblas-dev \
    #     libarpack2-dev \
    #     libparpack2-dev
    echo "Skipping system packages installation in favor of conda packages for building."
elif [ "${platform}" == "Darwin" ]; then
    echo "gfortran is expected to exist already."
fi

gfortran --version
which gfortran

# These packages are installed in the base environment but may be older
# versions. Explicitly upgrade them because they often create
# installation problems if out of date.
which python
python -VV

python -m pip install --upgrade pip setuptools wheel

# Generate .whl file.
if [ "${platform}" == "Linux" ]; then
    python setup.py sdist bdist_wheel --plat-name=manylinux2014_x86_64
else
    python setup.py sdist bdist_wheel
fi
ls -la dist/

# Install this package and the packages listed in requirements.txt.
python -m pip install dist/edrixs-*-cp${PYTHON_VERSION_NODOT}*.whl

# Install extra requirements for running tests and building docs.
python -m pip install -r requirements-dev.txt

# List the depencencies
python -m pip list
