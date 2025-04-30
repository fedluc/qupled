FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Install Python and clean up
RUN apt-get update && apt-get install -y \
    python3 python3-venv python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Set up Python environment
ENV PYTHON_VENV_ROOT=/workspaces/qupled_venv
RUN python3 -m venv $PYTHON_VENV_ROOT && \
    chown -R vscode:vscode $PYTHON_VENV_ROOT

# Copy project files, set permissions, and install dependencies
ENV SETUP_DIR=/setup_dir
COPY dev $SETUP_DIR/dev/
COPY devtool $SETUP_DIR/devtool
WORKDIR $SETUP_DIR
RUN chmod +x devtool && \
    . $PYTHON_VENV_ROOT/bin/activate && \
    ./devtool install-deps
WORKDIR /
RUN rm -rf $SETUP_DIR

# Automatically activate the Python environment on container start
RUN echo "source $PYTHON_VENV_ROOT/bin/activate" >> /etc/bash.bashrc