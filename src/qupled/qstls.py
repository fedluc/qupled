from __future__ import annotations
import numpy as np
from . import database
from . import native
from . import stls
from . import output


class Qstls(stls.Stls):
    """
    Class used to solve the Qstls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.Qstls
        self.native_inputs = native.QstlsInput()

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
        inputs = data[database.DataBaseHandler.INPUTS_TABLE_NAME]
        results = data[database.DataBaseHandler.RESULTS_TABLE_NAME]
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
        self.fixed: str = ""
        """ Name of the file storing the fixed component of the auxiliary density 
        response in the QSTLS scheme. """
        self.guess: Guess = Guess()
        """Initial guess. Default = ``Qstls.Guess()``"""
        # Undocumented default values
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
