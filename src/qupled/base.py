from __future__ import annotations

import os
import sys
import json
import sqlalchemy as sql
import sqlalchemy.orm as orm
import numpy as np
import io
import pandas as pd

from . import native
from . import util

DB_BASE = orm.declarative_base()


# -----------------------------------------------------------------------
# DataBaseHandler class
# -----------------------------------------------------------------------


class DataBaseHandler:

    TYPE_MAPPING = {int: sql.Integer, float: sql.Float, str: sql.String}

    def __init__(self, input_cls, input_table_name, result_cls, result_table_name):
        self.input_table = self._build_input_table(input_cls, input_table_name)
        self.result_table = self._build_result_table(result_cls, result_table_name)

    def _build_input_table(self, input_cls, input_table_name):
        columns = {
            "__tablename__": input_table_name,
            "id": sql.Column(sql.Integer, primary_key=True, autoincrement=True),
        }
        for attr, value in input_cls.__dict__.items():
            if not attr.startswith("__") and not callable(value):
                sql_type = self.TYPE_MAPPING.get(type(value), sql.String)
                columns[attr] = sql.Column(sql_type, nullable=False)
        return type("InputTable", (DB_BASE,), columns)

    def _build_result_table(self, result_cls, result_table_name):
        columns = {
            "__tablename__": result_table_name,
            "id": sql.Column(sql.Integer, primary_key=True, autoincrement=True),
        }
        for attr, value in result_cls.__dict__.items():
            if not attr.startswith("__") and not callable(value):
                columns[attr] = sql.Column(sql.LargeBinary, nullable=True)
        return type("ResultTable", (DB_BASE,), columns)


# -----------------------------------------------------------------------
# Input class
# -----------------------------------------------------------------------


class Input:

    @staticmethod
    def to_native(input_cls, native_input: any) -> any:
        for attr, value in input_cls.__dict__.items():
            setattr(native_input, attr, value)
        return native_input

    @staticmethod
    def to_database_table(input_cls, table_cls) -> any:
        input_data = {
            attr: (json.dumps(value) if isinstance(value, (list, dict)) else value)
            for attr, value in input_cls.__dict__.items()
        }
        return table_cls(**input_data)


# -----------------------------------------------------------------------
# Result class
# -----------------------------------------------------------------------


class Result:

    @staticmethod
    def from_native(result_cls, native_scheme: any) -> any:
        for attr in result_cls.__dict__.keys():
            value = getattr(native_scheme, attr)
            setattr(result_cls, attr, value) if value is not None else None

    @staticmethod
    def to_database_table(result_cls, table_cls) -> any:
        result_data = {
            attr: (Result.numpy_to_bytes(value) if value is not None else None)
            for attr, value in result_cls.__dict__.items()
        }
        return table_cls(**result_data)

    @staticmethod
    def numpy_to_bytes(arr: np.array) -> bytes:
        arr_bytes = io.BytesIO()
        np.save(arr_bytes, arr)
        return arr_bytes.getvalue()


# -----------------------------------------------------------------------
# ClassicScheme class
# -----------------------------------------------------------------------


class ClassicScheme:

    def __init__(self):
        # File to store output on disk
        self.hdf_file_name: str = None  #: Name of the output file.
        self.session = None

    # Compute the scheme
    def _compute(self, scheme) -> None:
        self.hdf_file_name = self._get_hdf_file(scheme.inputs)
        status = scheme.compute()
        self._check_status_and_clean(status, scheme.recovery)

    # Check that the dielectric scheme was solved without errors
    @util.MPI.run_only_on_root
    def _check_status_and_clean(self, status: bool, recovery: str) -> None:
        """Checks that the scheme was solved correctly and removes temporary files generated at run-time

        Args:
            status: status obtained from the native code. If status == 0 the scheme was solved correctly.
            recovery: name of the recovery file.
        """
        if status == 0:
            if os.path.isfile(recovery):
                os.remove(recovery)
            print("Dielectric theory solved successfully!")
        else:
            sys.exit("Error while solving the dielectric theory")

    # Save results to disk
    def _get_hdf_file(self, inputs) -> str:
        """Sets the name of the hdf file used to store the output

        Args:
            inputs: input parameters
        """
        coupling = inputs.coupling
        degeneracy = inputs.degeneracy
        theory = inputs.theory
        return f"rs{coupling:5.3f}_theta{degeneracy:5.3f}_{theory}.h5"

    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        inputs = scheme.inputs
        """Stores the results obtained by solving the scheme."""
        pd.DataFrame(
            {
                util.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                util.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                util.HDF.EntryKeys.THEORY.value: inputs.theory,
                util.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                util.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                util.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                util.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[util.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=util.HDF.EntryKeys.INFO.value, mode="w")
        if inputs.degeneracy > 0:
            pd.DataFrame(scheme.idr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.IDR.value
            )
            pd.DataFrame(scheme.sdr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.SDR.value
            )
            pd.DataFrame(scheme.slfc).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.SLFC.value
            )
        pd.DataFrame(scheme.ssf).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.SSF.value
        )
        pd.DataFrame(scheme.ssf_HF).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.SSF_HF.value
        )
        pd.DataFrame(scheme.wvg).to_hdf(
            self.hdf_file_name, key=util.HDF.EntryKeys.WVG.value
        )

    @util.MPI.run_only_on_root
    def _save_to_database(self, tables: list[any]) -> None:
        if self.session is None:
            engine = sql.create_engine("sqlite:///scheme_results.db")
            self.session = orm.sessionmaker(bind=engine)
        session = self.session()
        engine = session.get_bind()
        try:
            for table in tables:
                if not sql.inspect(engine).has_table(table.__tablename__):
                    table.__table__.create(engine)
                session.add(table)
                session.commit()
        except Exception as e:
            session.rollback()
            print(f"Database error: {e}")
        finally:
            session.close()

    # Compute radial distribution function
    def compute_rdf(
        self, rdf_grid: np.ndarray = None, write_to_hdf: bool = True
    ) -> np.array:
        """Computes the radial distribution function from the data stored in the output file.

        Args:
            rdf_grid: The grid used to compute the radial distribution function.
                Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`)
            write_to_hdf: Flag marking whether the rdf data should be added to the output hdf file, default to True

        Returns:
            The radial distribution function

        """
        if util.MPI().rank() > 0:
            write_to_hdf = False
        return util.HDF().compute_rdf(self.hdf_file_name, rdf_grid, write_to_hdf)

    # Compute the internal energy
    def compute_internal_energy(self) -> float:
        """Computes the internal energy from the data stored in the output file.

        Returns:
            The internal energy

        """
        return util.HDF().compute_internal_energy(self.hdf_file_name)

    # Plot results
    @util.MPI.run_only_on_root
    def plot(
        self,
        to_plot: list[str],
        matsubara: np.ndarray = None,
        rdf_grid: np.ndarray = None,
    ) -> None:
        """Plots the results stored in the output file.

        Args:
            to_plot: A list of quantities to plot. This list can include all the values written to the
                 output hdf file. The radial distribution function is computed and added to the output
                 file if necessary
            matsubara: A list of matsubara frequencies to plot. Applies only when the idr is plotted.
                (Default = None, all matsubara frequencies are plotted)
            rdf_grid: The grid used to compute the radial distribution function. Applies only when the radial
                distribution function is plotted. Default = ``None`` (see :func:`qupled.util.Hdf.computeRdf`).

        """
        if util.HDF.EntryKeys.RDF.value in to_plot:
            self.compute_rdf(rdf_grid)
        util.HDF().plot(self.hdf_file_name, to_plot, matsubara)


# -----------------------------------------------------------------------
# IterativeScheme class
# -----------------------------------------------------------------------


class IterativeScheme(ClassicScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> IterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = util.HDF().read(
            file_name, [util.HDF.EntryKeys.WVG.value, util.HDF.EntryKeys.SLFC.value]
        )
        return IterativeScheme.Guess(
            hdf_data[util.HDF.EntryKeys.WVG.value],
            hdf_data[util.HDF.EntryKeys.SLFC.value],
        )

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        inputs = scheme.inputs
        pd.DataFrame(
            {
                util.HDF.EntryKeys.COUPLING.value: inputs.coupling,
                util.HDF.EntryKeys.DEGENERACY.value: inputs.degeneracy,
                util.HDF.EntryKeys.ERROR.value: scheme.error,
                util.HDF.EntryKeys.THEORY.value: inputs.theory,
                util.HDF.EntryKeys.RESOLUTION.value: inputs.resolution,
                util.HDF.EntryKeys.CUTOFF.value: inputs.cutoff,
                util.HDF.EntryKeys.FREQUENCY_CUTOFF.value: inputs.frequency_cutoff,
                util.HDF.EntryKeys.MATSUBARA.value: inputs.matsubara,
            },
            index=[util.HDF.EntryKeys.INFO.value],
        ).to_hdf(self.hdf_file_name, key=util.HDF.EntryKeys.INFO.value)

    # Initial guess
    class Guess:

        def __init__(self, wvg: np.ndarray = None, slfc: np.ndarray = None):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.slfc = slfc
            """ Static local field correction. Default = ``None``"""

        def to_native(self) -> native.StlsGuess:
            native_guess = native.StlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess


# -----------------------------------------------------------------------
# QuantumIterativeScheme class
# -----------------------------------------------------------------------


class QuantumIterativeScheme(IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(file_name: str) -> QuantumIterativeScheme.Guess:
        """Constructs an initial guess object by extracting the information from an output file.

        Args:
            file_name : name of the file used to extract the information for the initial guess.
        """
        hdf_data = util.HDF().read(
            file_name,
            [
                util.HDF.EntryKeys.WVG.value,
                util.HDF.EntryKeys.SSF.value,
                util.HDF.EntryKeys.ADR.value,
                util.HDF.EntryKeys.MATSUBARA.value,
            ],
        )
        return QuantumIterativeScheme.Guess(
            hdf_data[util.HDF.EntryKeys.WVG.value],
            hdf_data[util.HDF.EntryKeys.SSF.value],
            np.ascontiguousarray(hdf_data[util.HDF.EntryKeys.ADR.value]),
            hdf_data[util.HDF.EntryKeys.MATSUBARA.value],
        )

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme)
        if scheme.inputs.degeneracy > 0:
            pd.DataFrame(scheme.adr).to_hdf(
                self.hdf_file_name, key=util.HDF.EntryKeys.ADR.value
            )

    # Initial guess
    class Guess:

        def __init__(
            self,
            wvg: np.ndarray = None,
            ssf: np.ndarray = None,
            adr: np.ndarray = None,
            matsubara: int = 0,
        ):
            self.wvg = wvg
            """ Wave-vector grid. Default = ``None``"""
            self.ssf = ssf
            """ Static structure factor. Default = ``None``"""
            self.adr = adr
            """ Auxiliary density response. Default = ``None``"""
            self.matsubara = matsubara
            """ Number of matsubara frequencies. Default = ``0``"""

        def to_native(self) -> native.QStlsGuess:
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess
