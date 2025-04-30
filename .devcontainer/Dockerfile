FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Install Python and clean up
RUN apt-get update && apt-get install -y \
    python3 python3-venv python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Set up Python environment
ENV PYTHON_VENV_ROOT=/workspaces/venv
RUN python3 -m venv $PYTHON_VENV_ROOT && \
    . $PYTHON_VENV_ROOT/bin/activate

# Copy project files and set permissions
COPY dev /workspaces/devtool/dev/
COPY devtool /workspaces/devtool/devtool
WORKDIR /workspaces/devtool
RUN chmod +x devtool
RUN . $PYTHON_VENV_ROOT/bin/activate && ./devtool install-deps