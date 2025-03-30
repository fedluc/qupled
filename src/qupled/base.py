from __future__ import annotations

import numpy as np

from . import native
from . import util
from qupled.database import DataBaseHandler

# -----------------------------------------------------------------------
# QuantumIterativeScheme class
# -----------------------------------------------------------------------


class QuantumIterativeScheme(IterativeScheme):

    # Set the initial guess from a dataframe produced in output
    @staticmethod
    def get_initial_guess(run_id: str, database_name: str | None = None) -> Guess:
        """
        Generates an initial guess for a quantum iterative scheme based on data
        retrieved from a database.

        Args:
            run_id (str): The unique identifier for the run whose data is to be retrieved.
            database_name (str | None, optional): The name of the database to query.
                If None, a default database is used.

        Returns:
            QuantumIterativeScheme.Guess: An object containing the initial guess values
            for the quantum iterative scheme, including results and inputs.

        Raises:
            Exception: If there is an issue reading data from the database.
        """
        result_names = ["wvg", "ssf", "adr"]
        input_names = ["matsubara"]
        data = util.DataBase.read_run(run_id, database_name, input_names, result_names)
        inputs = data[DataBaseHandler.INPUTS_TABLE_NAME]
        results = data[DataBaseHandler.RESULTS_TABLE_NAME]
        return QuantumIterativeScheme.Guess(
            results[result_names[0]],
            results[result_names[1]],
            results[result_names[2]],
            inputs[input_names[0]],
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
            """
            Converts the current object to a native `QStlsGuess` object.

            This method creates an instance of `native.QStlsGuess` and populates its
            attributes with the corresponding values from the current object's
            attributes. If an attribute's value is `None`, it is replaced with an
            empty NumPy array.

            Returns:
                native.QStlsGuess: A new instance of `native.QStlsGuess` with attributes
                populated from the current object.
            """
            native_guess = native.QstlsGuess()
            for attr, value in self.__dict__.items():
                native_value = value if value is not None else np.empty(0)
                setattr(native_guess, attr, native_value)
            return native_guess
