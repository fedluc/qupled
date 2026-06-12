from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

# Temporary files used for MPI executions
INPUT_FILE = Path("input.json")
RESULT_FILE = Path("results.json")
STATUS_FILE = Path("status.json")


class WorkerFiles:
    """Files used for MPI worker input, output, and status exchange.

    Args:
        directory: Optional directory used to hold the worker I/O files. Defaults
            to a temporary directory owned by these worker files.
    """

    def __init__(self, directory: Path | None = None):
        self._owns_directory = directory is None
        if directory is None:
            directory = Path(tempfile.mkdtemp(prefix="qupled-mpi-"))
        self.directory = directory
        self.input_file = self.directory / INPUT_FILE.name
        self.result_file = self.directory / RESULT_FILE.name
        self.status_file = self.directory / STATUS_FILE.name

    def cleanup(self):
        """Clean up the owned temporary worker directory."""
        if self._owns_directory and self.directory.exists():
            shutil.rmtree(self.directory)

    def _write_json_atomically(self, target_file: Path, data):
        """Write JSON data through a temporary file and atomic replace.

        Args:
            target_file: Final JSON file path.
            data: JSON-serializable data to write.
        """
        temp_file = None
        try:
            with tempfile.NamedTemporaryFile(
                "w",
                delete=False,
                dir=target_file.parent,
                prefix=f".{target_file.name}.",
                suffix=".tmp",
            ) as f:
                temp_file = Path(f.name)
                json.dump(data, f)
            temp_file.replace(target_file)
        except Exception:
            if temp_file is not None and temp_file.exists():
                temp_file.unlink()
            raise

    def write_inputs(self, inputs):
        """Write input data to the MPI worker input file.

        Args:
            inputs: Serializable input object for the solver run.
        """
        with self.input_file.open("w") as f:
            json.dump(inputs.to_dict(), f)

    def read_inputs(self, input_cls):
        """Read input data from the MPI worker input file.

        Args:
            input_cls: Serializable input class used to reconstruct the inputs.

        Returns:
            Reconstructed input object.
        """
        with self.input_file.open() as f:
            input_dict = json.load(f)
        return input_cls.from_dict(input_dict)

    def write_results(self, scheme, result_cls):
        """Write solver results from the root MPI rank.

        Args:
            scheme: Native scheme instance.
            result_cls: Serializable result class used to convert native results.
        """
        if scheme.is_root:
            results = result_cls()
            results.from_native(scheme)
            self._write_json_atomically(self.result_file, results.to_dict())

    def read_results(self, result_cls):
        """Read solver results from the MPI worker result file.

        Args:
            result_cls: Serializable result class used to reconstruct the results.

        Returns:
            Reconstructed result object.
        """
        with self.result_file.open() as f:
            result_dict = json.load(f)
        return result_cls.from_dict(result_dict)

    def write_status(self, scheme, status):
        """Write solver status from the root MPI rank.

        Args:
            scheme: Native scheme instance.
            status: Solver status to serialize.
        """
        if scheme.is_root:
            self._write_json_atomically(self.status_file, status)

    def read_status(self):
        """Read solver status from the MPI worker status file.

        Returns:
            Deserialized solver status.
        """
        with self.status_file.open() as f:
            status = json.load(f)
        return status
