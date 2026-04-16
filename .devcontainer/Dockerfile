FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Install Python and clean up
RUN apt-get update && apt-get install -y \
    curl python3 python3-venv python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install uv and configure the project environment location
RUN curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh
ENV PYTHON_VENV_ROOT=/workspaces/qupled_venv
ENV UV_PROJECT_ENVIRONMENT=$PYTHON_VENV_ROOT

# Copy project files, set permissions, and install dependencies
ENV SETUP_DIR=/setup_dir
ENV DEVTOOL_SKIP_PROJECT_INSTALL=1
COPY dev $SETUP_DIR/dev/
COPY pyproject.toml $SETUP_DIR/pyproject.toml
COPY devtool $SETUP_DIR/devtool
WORKDIR $SETUP_DIR
RUN chmod +x devtool && \
    ./devtool install-deps && \
    chown -R vscode:vscode $PYTHON_VENV_ROOT
WORKDIR /
RUN rm -rf $SETUP_DIR
