FROM ubuntu:24.04 AS base

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

FROM base AS dep_builder

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

ARG ARPACK_VERSION=3.9.1
ARG ARPACK_SHA256=f6641deb07fa69165b7815de9008af3ea47eb39b2bb97521fbf74c97aba6e844
RUN curl --fail --show-error --location --retry 3 \
      "https://github.com/opencollab/arpack-ng/archive/refs/tags/${ARPACK_VERSION}.tar.gz" \
      --output /tmp/arpack-ng.tar.gz && \
  echo "${ARPACK_SHA256}  /tmp/arpack-ng.tar.gz" | sha256sum --check - && \
  tar --extract --gzip --file /tmp/arpack-ng.tar.gz --directory /tmp && \
  cmake -S "/tmp/arpack-ng-${ARPACK_VERSION}" \
      -B /tmp/arpack-ng-build \
      -DMPI=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/edrixs && \
  cmake --build /tmp/arpack-ng-build --parallel 4 && \
  ctest --test-dir /tmp/arpack-ng-build --output-on-failure && \
  cmake --install /tmp/arpack-ng-build && \
  rm -rf /tmp/arpack-ng.tar.gz "/tmp/arpack-ng-${ARPACK_VERSION}" \
      /tmp/arpack-ng-build

FROM dep_builder AS app_builder

COPY --chown=docker docker/requirements.txt /tmp/edrixs-docker-requirements.txt
RUN python -m pip install --no-cache-dir --upgrade \
      pip==26.2.1 setuptools==84.0.0 wheel==0.48.0 && \
  python -m pip install --no-cache-dir \
      --requirement /tmp/edrixs-docker-requirements.txt

COPY --chown=docker . edrixs

RUN export VERBOSE=1 FFLAGS="-Wall -Wunused -Wextra -Wno-maybe-uninitialized -Ofast -faggressive-loop-optimizations -fno-tree-pre" && \
  python -m pip install --no-build-isolation --no-deps --verbose ./edrixs

FROM app_builder AS developer
RUN python -m pip install --no-cache-dir \
    --requirement edrixs/requirements-dev.txt

FROM base AS app

# see https://github.com/open-mpi/ompi/issues/4948
ENV LD_LIBRARY_PATH="/opt/edrixs/lib" \
    CMAKE_PREFIX_PATH="/opt/edrixs" \
    OMPI_MCA_btl_vader_single_copy_mechanism="none"

COPY --from=app_builder /opt /opt
COPY --from=app_builder --chown=docker /home/docker/edrixs/examples examples
