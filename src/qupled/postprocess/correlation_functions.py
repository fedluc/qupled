import numpy as np

from qupled.database.base_tables import ConflictMode
from qupled.database.database_handler import DataBaseHandler
from qupled.postprocess.output import DataBase, OutputType
from qupled.schemes import hf, rpa
from qupled.util.dimension import Dimension


class CorrelationFunctions:
    """Post-processing module to compute correlation functions from existing scheme runs."""

    def __init__(self, database_name: str | None = None):
        self.database_name = database_name

    def compute_rdf(self, run_id: int, rdf_grid: np.ndarray | None = None) -> None:
        """
        Compute the radial distribution function (RDF) for an existing scheme run
        and save it back to the database.

        Args:
            run_id: The ID of the scheme run.
            rdf_grid: Grid points at which to evaluate the RDF.
                If not provided, defaults to np.arange(0.0, 10.0, 0.01).
        """
        data = DataBase.read_run(
            run_id,
            type=OutputType.SCHEME,
            database_name=self.database_name,
            input_names=["dimension"],
            result_names=["wvg", "ssf"],
        )
        dimension = Dimension.from_dict(data.inputs["dimension"])
        result = hf.Result.from_dict(data.results)
        result.compute_rdf(dimension, rdf_grid)
        db_handler = DataBaseHandler(self.database_name)
        db_handler.scheme_tables.run_id = run_id
        db_handler.scheme_tables.insert_results(
            {"rdf": result.rdf, "rdf_grid": result.rdf_grid},
            conflict_mode=ConflictMode.UPDATE,
        )

    def compute_itcf(self, run_id: int, tau: np.ndarray | None = None) -> None:
        """
        Compute the imaginary-time correlation function (ITCF) for an existing scheme run
        and save it back to the database.

        For HF runs, uses the non-interacting formula; all other theories use the
        interacting formula.

        The imaginary-time grid is in absolute units [0, 1/theta]. If tau is not provided,
        a default grid is generated based on the degeneracy parameter of the run (see
        :meth:`qupled.schemes.hf.Result.compute_itcf` for details). For finite-temperature
        runs, tau values exceeding 1/theta are silently discarded.

        Args:
            run_id: The ID of the scheme run.
            tau: Imaginary-time grid points in absolute units [0, 1/theta] at which to
                evaluate the ITCF. If not provided, a default grid is used.
        """
        data = DataBase.read_run(
            run_id,
            type=OutputType.SCHEME,
            database_name=self.database_name,
            result_names=["wvg", "idr", "chemical_potential", "lfc"],
        )
        inputs = hf.Input.from_dict(data.inputs)
        if data.run["theory"] == "HF":
            results = hf.Result.from_dict(data.results)
        else:
            results = rpa.Result.from_dict(data.results)
        results.compute_itcf(inputs, tau)
        db_handler = DataBaseHandler(self.database_name)
        db_handler.scheme_tables.run_id = run_id
        db_handler.scheme_tables.insert_results(
            {"itcf": results.itcf, "tau": results.tau},
            conflict_mode=ConflictMode.UPDATE,
        )
