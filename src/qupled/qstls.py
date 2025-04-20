from __future__ import annotations

import numpy as np

from . import database
from . import native
from . import output
from . import stls


class Qstls(stls.Stls):
    """
    Class used to solve the Qstls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.Qstls
        self.native_inputs_cls = native.QstlsInput

    def compute(self, inputs: Input):
        self.find_fixed_adr_in_database(inputs)
        super().compute(inputs)

    def find_fixed_adr_in_database(self, inputs: Input):
        """
        Searches the database for a run with matching parameters and assigns its ID to the input object.

        This method iterates through all runs in the database and checks if a run matches the given
        input parameters (degeneracy, theory, cutoff, matsubara, and resolution). If a match is found,
        the `fixed_run_id` attribute of the input object is updated with the corresponding run ID.

        Args:
            inputs (Input): The input parameters.

        Returns:
            None: The method updates the `fixed_run_id` attribute of the `inputs` object if a match is found.
        """
        runs = self.db_handler.inspect_runs()
        inputs.fixed_run_id = None
        for run in runs:
            database_keys = database.DataBaseHandler.TableKeys
            same_degeneracy = run[database_keys.DEGENERACY.value] == inputs.degeneracy
            same_theory = run[database_keys.THEORY.value] == inputs.theory
            if not same_theory or not same_degeneracy:
                continue
            run_id = run[database_keys.PRIMARY_KEY.value]
            run_inputs = self.db_handler.get_inputs(run_id)
            if (
                run_inputs["cutoff"] == inputs.cutoff
                and run_inputs["matsubara"] == inputs.matsubara
                and run_inputs["resolution"] == inputs.resolution
            ):
                print(f"Loading fixed ADR from database for run_id = {run_id}")
                inputs.fixed_run_id = run_id
                return

    @staticmethod
    def get_initial_guess(run_id: str, database_name: str | None = None) -> Guess:
        """
        Retrieves the initial guess for a computation based on a specific run ID
        from a database.

        Args:
            run_id: The unique identifier for the run whose data is to be retrieved.
            database_name: The name of the database to query.
                If None, the default database is used.

        Returns:
            Guess: An object containing the initial guess values, including results
            and inputs extracted from the database.
        """
        result_names = ["wvg", "ssf", "adr"]
        input_names = ["matsubara"]
        data = output.DataBase.read_run(
            run_id, database_name, input_names, result_names
        )
        inputs = data[database.DataBaseHandler.INPUT_TABLE_NAME]
        results = data[database.DataBaseHandler.RESULT_TABLE_NAME]
        return Guess(
            results[result_names[0]],
            results[result_names[1]],
            results[result_names[2]],
            inputs[input_names[0]],
        )


# Input class
class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qstls.Qstls` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        super().__init__(coupling, degeneracy)
        self.guess: Guess = Guess()
        """Initial guess. Default = ``qstls.Guess()``"""
        # Undocumented default values
        self.fixed_run_id: int | None = None
        self.theory = "QSTLS"


class Result(stls.Result):
    """
    Class used to store the results for the :obj:`qupled.qstls.Qstls` class.
    """

    def __init__(self):
        super().__init__()
        self.adr = None
        """Auxiliary density response"""


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
            if value is not None:
                setattr(native_guess, attr, value)
        return native_guess
