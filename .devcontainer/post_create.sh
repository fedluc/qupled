#!/usr/bin/env bash
set -euo pipefail

# Ensure the user-local bin directory is available in this shell.
export PATH="${HOME}/.local/bin:${PATH}"

# Install uv if the base image does not already provide it.
if ! command -v uv >/dev/null 2>&1; then
  mkdir -p "${HOME}/.local/bin"
  curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="${HOME}/.local/bin" sh
fi

# Create the repository-local virtual environment used by the devcontainer.
uv venv .venv --allow-existing

# Install foga into that environment so the bootstrap commands below can use it.
uv pip install --python .venv/bin/python foga

# Activate the environment for the remaining setup steps.
. .venv/bin/activate

# Install required system packages for native builds and docs generation.
foga install --target apt-get
foga install --target dev-python
