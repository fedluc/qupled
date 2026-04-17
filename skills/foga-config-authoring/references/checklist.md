# Foga Config Checklist

## Inspect First

Read only the files that define or constrain workflow behavior:

- `pyproject.toml`
- `setup.py`
- `tox.ini`
- `.github/workflows/*`
- `CMakeLists.txt` and included CMake modules
- build or maintenance scripts under `dev/`, `scripts/`, `tools/`, or similar
- test layout under `tests/`
- docs or README sections that document installation, build, or test commands

## Build Mapping

### Python Build

Use `build.python.backend: python-build` when the package is built with `python -m build`.

Add:

- `env` when packaging relies on environment variables
- `args` only for extra `python -m build` flags

### C++ Build

Use `build.cpp.backend: cmake` when the repository has an explicit standalone CMake build.

Set:

- `source_dir`
- `build_dir`
- optional `generator`
- `configure_args`
- `build_args`
- `targets` only when a stable named target should be the default
- `build_args` when resource limits or tool behavior require explicit parallelism or other build flags

Do not collapse standalone C++ builds into `python-build` just because the package also compiles a compiled extension during wheel creation. If the repo documents stable top-level CMake targets, prefer those targets over inventing a `ctest` runner.

## Test Mapping

### Pytest

Use a `pytest` runner when the real command is fundamentally pytest.

Set:

- `path`
- `marker` when the repository uses marker-separated suites
- `args` for verbosity or other flags

Split runners by behavior when the repo has distinct suites such as `unit`, `cpp`, or `integration`.

### Tox

Use a `tox` runner only when tox is materially part of the workflow:

- environment matrix behavior matters
- install isolation is intentional
- dependencies or commands differ in tox-specific ways
- reproducing the behavior directly in `foga` would be less clear

If tox only shells out to pytest with minor setup, usually prefer `pytest` plus hooks.

### CTest

Use `ctest` for standalone C++ test binaries.

Set:

- `build_dir` always
- `source_dir` when `foga` must configure/build before `ctest`
- `configure_args` for feature toggles needed by tests
- `target` when the repo builds a dedicated test target before `ctest`
- `args` for `-V` or other ctest flags

## Hooks

Use hooks sparingly.

Good uses:

- stage example files needed by integration tests
- create or remove temporary directories the test suite assumes
- reproduce one small tox or CI preparation step
- bridge one missing repo-specific setup step that is not already handled by the devcontainer or project bootstrap

Bad uses:

- wrapping the entire build or test workflow in custom shell commands
- hiding the real backend command from `foga`
- repeatedly reinstalling dependencies in hooks when the environment is already provisioned

## Profiles

Add profiles only for real modes present in the target repo, such as:

- MPI on/off
- release vs debug C++ options
- platform-specific environment variables
- alternate deploy targets

Prefer profile overrides when the same repository supports multiple valid workflows and the difference is not just one ad hoc command invocation. Keep example configs small: add only the profiles that correspond to documented or stable variants.

## Final Checks

Before finishing:

1. Confirm every configured backend exists in `foga`.
2. Confirm each path is relative to the target repository root where `foga.yml` will live.
3. Validate with `foga validate`.
4. Inspect with `foga inspect`.
5. Explain any assumptions that could still depend on the target machine, such as external libraries, platform-specific env vars, or memory-sensitive build parallelism.
6. If profiles were added, inspect at least one representative profile to confirm the override merges as intended.
7. Keep `clean.paths` aligned with any profile-specific build directories you introduced.
