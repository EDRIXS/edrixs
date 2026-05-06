FROM ubuntu:24.04 as base

# general environment for docker
ENV DEBIAN_FRONTEND=noninteractive

# base packages
RUN apt-get update && apt-get install -y --no-install-recommends sudo curl vim libgfortran-13-dev libopenmpi-dev libopenblas-dev python3 libpython3-dev python3-pip python3-venv && rm -rf /var/lib/apt/lists/*

# create docker user
RUN useradd -m -s /bin/bash -u 1999 docker && echo "docker:docker" | chpasswd && adduser docker sudo

# enable passwordless sudo
RUN echo "docker ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/docker

# allow docker user to install into /opt
RUN sudo chown docker:docker /opt

USER    docker
WORKDIR /home/docker

# create venv
ENV VIRTUAL_ENV=/opt/edrixs
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

FROM base as dep_builder

ENV LD_LIBRARY_PATH="/opt/edrixs/lib" \
    CMAKE_PREFIX_PATH="/opt/edrixs"

# install build tools
RUN sudo apt-get update && sudo apt-get install -y --no-install-recommends \
    autoconf \
    build-essential \
    gfortran \
    ca-certificates \
    coreutils \
    curl \
    git \
    cmake

RUN curl -L https://github.com/opencollab/arpack-ng/archive/refs/tags/3.8.0.tar.gz | tar xvz && \
  cd arpack-ng-3.8.0 && \
  mkdir build && \
  cd build && \
  cmake -DMPI=ON -DCMAKE_INSTALL_PREFIX=/opt/edrixs .. && \
  make -j4 && \
  make test && \
  make install

FROM dep_builder as app_builder

RUN pip install --upgrade pip setuptools
RUN pip install numpy scipy sympy matplotlib sphinx mpi4py ipython jupyter jupyterlab ipympl ipywidgets lmfit

COPY --chown=docker . edrixs

RUN export VERBOSE=1 FFLAGS="-Wall -Wunused -Wextra -Wno-maybe-uninitialized -Ofast -faggressive-loop-optimizations -fno-tree-pre" && \
  pip install -v ./edrixs

FROM app_builder as developer
  RUN pip install -r edrixs/requirements-dev.txt

FROM base as app

# see https://github.com/open-mpi/ompi/issues/4948
ENV LD_LIBRARY_PATH="/opt/edrixs/lib" \
    CMAKE_PREFIX_PATH="/opt/edrixs" \
    OMPI_MCA_btl_vader_single_copy_mechanism="none"

COPY --from=app_builder /opt /opt
COPY --from=app_builder --chown=docker /home/docker/edrixs/examples examples
