import numpy as np

from qupled.database.base_tables import ConflictMode
from qupled.database.database_handler import DataBaseHandler
from qupled.postprocess.output import DataBase, OutputType
from qupled.schemes import hf
from qupled.util.dimension import Dimension


def compute_rdf(
    run_id: int, rdf_grid: np.ndarray | None = None, database_name: str | None = None
) -> None:
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
        database_name=database_name,
        input_names=["dimension"],
        result_names=["wvg", "ssf"],
    )
    dimension = Dimension.from_dict(data.inputs["dimension"])
    result = hf.Result.from_dict(data.results)
    result.compute_rdf(dimension, rdf_grid)
    db_handler = DataBaseHandler(database_name)
    db_handler.scheme_tables.run_id = run_id
    db_handler.scheme_tables.insert_results(
        {"rdf": result.rdf, "rdf_grid": result.rdf_grid},
        conflict_mode=ConflictMode.UPDATE,
    )


def compute_itcf(
    run_id: int, tau: np.ndarray | None = None, database_name: str | None = None
) -> None:
    """
    Compute the imaginary-time correlation function (ITCF) for an existing scheme run
    and save it back to the database.

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
        database_name=database_name,
        result_names=["wvg", "idr", "chemical_potential", "lfc"],
    )
    inputs = hf.Input.from_dict(data.inputs)
    results = hf.Result.from_dict(data.results)
    results.compute_itcf(inputs, tau)
    db_handler = DataBaseHandler(database_name)
    db_handler.scheme_tables.run_id = run_id
    db_handler.scheme_tables.insert_results(
        {"itcf": results.itcf, "tau": results.tau},
        conflict_mode=ConflictMode.UPDATE,
    )


def compute_dsf(
    run_id: int,
    frequency: np.ndarray | None = None,
    database_name: str | None = None,
) -> None:
    """
    Compute the dynamic structure factor for an existing static-scheme run and
    save it back to the database.

    Args:
        run_id: The ID of the scheme run.
        frequency: Non-negative real-frequency grid. If not provided, a default
            grid is generated from the run inputs.
        database_name: Name of the database to read from.
    """
    data = DataBase.read_run(
        run_id,
        type=OutputType.SCHEME,
        database_name=database_name,
        result_names=["wvg", "chemical_potential", "lfc"],
    )
    inputs = hf.Input.from_dict(data.inputs)
    results = hf.Result.from_dict(data.results)
    results.compute_dsf(inputs, frequency)
    db_handler = DataBaseHandler(database_name)
    db_handler.scheme_tables.run_id = run_id
    db_handler.scheme_tables.insert_results(
        {"dsf": results.dsf, "frequency": results.frequency},
        conflict_mode=ConflictMode.UPDATE,
    )


def compute_itcf_from_dsf(
    run_id: int, tau: np.ndarray | None = None, database_name: str | None = None
) -> None:
    """
    Compute the imaginary-time correlation function from a stored dynamic
    structure factor and save it back to the database.

    Args:
        run_id: The ID of the scheme run.
        tau: Imaginary-time grid in absolute units [0, 1/theta].
        database_name: Name of the database to read from.
    """
    data = DataBase.read_run(
        run_id,
        type=OutputType.SCHEME,
        database_name=database_name,
        result_names=["frequency", "dsf", "ssf"],
    )
    inputs = hf.Input.from_dict(data.inputs)
    results = hf.Result.from_dict(data.results)
    results.compute_itcf_from_dsf(inputs, tau)
    db_handler = DataBaseHandler(database_name)
    db_handler.scheme_tables.run_id = run_id
    db_handler.scheme_tables.insert_results(
        {
            "itcf_from_dsf": results.itcf_from_dsf,
            "tau": results.tau,
            "ssf_from_dsf": results.ssf_from_dsf,
            "delta_ssf_from_dsf": results.delta_ssf_from_dsf,
        },
        conflict_mode=ConflictMode.UPDATE,
    )
