---
name: foga-config-authoring
description: Author or update `foga.yml` files for repositories that use Python, C++ build systems, test runners, packaging, and deployment workflows. Use when Codex should inspect another repository, map its real build and test commands onto `foga` backends such as `cmake`, `python-build`, `pytest`, `tox`, `ctest`, and `twine`, and produce a validated, reusable `foga.yml` instead of guessing.
---

# Foga Config Authoring

Use this skill to translate an existing repository workflow into a working `foga.yml`.

## Workflow

1. Confirm the available `foga` schema and backends from the environment you are running in.
   - First prefer local `foga` help, examples, installed package files, or repository source if they are available.
   - If working inside the `foga` repository, reading `README.md`, `src/foga/config.py`, and `examples/*/foga.yml` is appropriate.
   - If working outside the `foga` repository, inspect the installed package, local docs, or `foga validate` and `foga inspect --help` output instead of assuming repo-local source files exist.
   - The published documentation is available at https://fedluc.github.io/foga/.
2. Inspect the target repository before writing config. Start with `pyproject.toml`, `setup.py`, `tox.ini`, build scripts, CI workflows, `CMakeLists.txt`, and tests.
3. Identify the real workflows separately:
   - Python package build
   - C++ build
   - Python tests
   - C++ tests
   - Deploy steps, if requested
4. Map each workflow to the narrowest `foga` backend that expresses it directly. Prefer `pytest`, `ctest`, `cmake`, and `python-build` over wrapper layers when those wrappers add little value.
5. Add hooks only for repository-specific setup that the direct backend cannot express, such as copying fixture files or preparing working directories. Prefer zero-hook configs when the environment is already provisioned by a devcontainer, CI image, or repo bootstrap step.
6. Add profiles only when the target repo has real environment variants such as MPI, debug/release toggles, platform-specific environment variables, alternate C++ options, or documented named target subsets. Keep the profile set minimal for example configs.
7. Validate the generated file with the locally available `foga` CLI before finishing.

## Authoring Rules

- Do not infer commands from project type alone. Derive them from the target repo's actual files.
- Keep Python and C++ workflows separate when the repository treats them separately.
- Use `build.cpp.targets` only when the repository has a stable named target worth selecting explicitly. If the project documents top-level CMake targets such as `check`, `cpptest`, `pytest`, or `test_cmake_build`, prefer modelling those targets directly instead of inventing a `ctest` runner.
- For `ctest`, set `source_dir`, `build_dir`, configure arguments, and optional `target` when test binaries must be prepared before `ctest` runs.
- Prefer `pytest` runners over `tox` when `tox` mostly forwards to pytest and `foga` can represent the same behavior directly.
- Keep environment variables close to the workflow that needs them. Use profiles when the values vary by mode rather than by command type.
- If tests rely on side effects from CI or tox setup, reproduce only the minimal missing behavior with hooks. If the repo already provisions its environment separately, remove setup hooks rather than encoding package installation into every workflow.

## Validation

- Prefer the installed CLI entrypoint when available, for example `foga validate --config <path>/foga.yml` or the equivalent local invocation supported by the environment.
- Run a validation command against the generated file.
- Run an inspect command against the generated file.
- If possible, note the intended `foga build ...` and `foga test ...` commands for the target repo.

## References

- Read [references/checklist.md](references/checklist.md) for the inspection checklist and mapping heuristics.

## Work Log

- 2026-04-13: Derived `pybind11/foga.yml` from the checked-in contributing guide, `CMakePresets.json`, and test targets. Recorded three reusable heuristics: prefer direct CMake targets over `ctest` when the repo exposes stable targets like `check` or `cpptest`; start with a hook-free config when the environment is already provisioned; and keep profiles minimal but aligned with documented variants such as `release`, `venv`, or `tidy`.

- 2026-04-12: Derived `examples/qupled/foga.yml` by inspecting `qupled` build scripts, `tox.ini`, CMake files, and pytest markers instead of assuming a generic Python package layout.
- 2026-04-12: Recorded reusable guidance from that work: prefer direct backends over thin wrappers, model standalone C++ test builds separately from extension-module builds when the repo does, and use hooks only for missing setup behavior such as integration-test fixture staging.
