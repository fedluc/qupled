FROM quay.io/pypa/manylinux_2_28_x86_64

# Install dev tools and system dependencies
RUN dnf -y update && \
    dnf -y groupinstall "Development Tools" && \
    dnf -y install \
        openmpi \
        openmpi-devel \
        gsl-devel \
        sqlite-devel \
        git \
        wget \
        cmake \
    && dnf clean all

# Initialize environment-modules in login shells
RUN echo 'module load mpi/openmpi-x86_64' >> /etc/bashrc

# Set Python to 3.10 (or another available version under /opt/python)
ENV PYTHON_ROOT=/opt/python/cp310-cp310
ENV PATH="${PYTHON_ROOT}/bin:$PATH"
ENV LD_LIBRARY_PATH="${PYTHON_ROOT}/lib:$LD_LIBRARY_PATH"

# Build SQLiteCpp from source
RUN cd /tmp && \
    git clone https://github.com/SRombauts/SQLiteCpp.git && \
    cd SQLiteCpp && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd / && rm -rf /tmp/SQLiteCpp

# Set working directory and copy source
WORKDIR /src

# Set entrypoint to ensure MPI module is loaded
SHELL ["/bin/bash", "-lc"]
