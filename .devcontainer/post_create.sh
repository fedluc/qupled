#!/usr/bin/env bash
set -euo pipefail

# Ensure the user-local bin directory is available in this shell.
export PATH="${HOME}/.local/bin:${PATH}"

# Install uv if the base image does not already provide it.
if ! command -v uv >/dev/null 2>&1; then
  mkdir -p "${HOME}/.local/bin"
  curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="${HOME}/.local/bin" sh
fi

# Sync the repository-local environment used by the devcontainer.
uv sync --group dev --no-install-project

# Install required system packages for native builds and docs generation.
. .venv/bin/activate
foga install --target apt-get
foga install --target dev-env
