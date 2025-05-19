from __future__ import annotations

import numpy as np

from . import hf
from . import native
from . import output
from . import rpa


class Stls(rpa.Rpa):
    """
    Class used to solve the Stls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.Stls
        self.native_inputs_cls = native.StlsInput

    @staticmethod
    def get_initial_guess(run_id: int, database_name: str | None = None) -> Guess:
        """Constructs an initial guess object by extracting the information from a database.

        Args:
            run_id: The unique identifier for the run whose data is to be retrieved.
            database_name: The name of the database to query.
                If None, the default database will be used.

        Returns:
            An instance of Guess containing the initial guess data.
        """
        names = ["wvg", "ssf"]
        data = output.DataBase.read_results(run_id, database_name, names)
        return Guess(data[names[0]], data[names[1]])


class Input(rpa.Input):
    """
    Class used to manage the input for the :obj:`qupled.stls.Stls` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        super().__init__(coupling, degeneracy)
        self.error: float = 1.0e-5
        """Minimum error for convergence. Default = ``1.0e-5``"""
        self.mixing: float = 1.0
        """Mixing parameter. Default = ``1.0``"""
        self.iterations: int = 1000
        """Maximum number of iterations. Default = ``1000``"""
        self.guess: Guess = Guess()
        """Initial guess. Default = ``stls.Guess()``"""
        # Undocumented default values
        self.theory: str = "STLS"


class Result(hf.Result):
    """
    Class used to store the results for the :obj:`qupled.stls.Stls` class.
    """

    def __init__(self):
        super().__init__()
        self.error: float = None
        """Residual error in the solution"""


class Guess:

    def __init__(self, wvg: np.ndarray = None, ssf: np.ndarray = None):
        self.wvg = wvg
        """ Wave-vector grid. Default = ``None``"""
        self.ssf = ssf
        """ Static structure factor. Default = ``None``"""

    def to_native(self) -> native.Guess:
        """
        Converts the current object to a native `Guess` object.

        This method iterates over the attributes of the current object and
        assigns their values to a new `Guess` object. If an attribute's
        value is `None`, it is replaced with an empty NumPy array.

        Returns:
            native.StlsGuess: A new instance of `Guess` with attributes
            copied from the current object.
        """
        native_guess = native.Guess()
        for attr, value in self.__dict__.items():
            if value is not None:
                setattr(native_guess, attr, value)
        return native_guess
