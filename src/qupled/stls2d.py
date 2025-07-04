from __future__ import annotations

from dataclasses import field
import numpy as np

from . import hf2d
from . import native
from . import output
from . import rpa2d
from . import serialize


class Solver(rpa2d.Solver):
    """
    Class used to solve the Stls scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.Stls
    native_inputs_cls = native.StlsInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()

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


@serialize.serializable_dataclass
class Input(rpa2d.Input):
    """
    Class used to manage the input for the :obj:`qupled.stls.Stls` class.
    """

    error: float = 1.0e-5
    """Minimum error for convergence. Default = ``1.0e-5``"""
    mixing: float = 1.0
    """Mixing parameter. Default = ``1.0``"""
    iterations: int = 1000
    """Maximum number of iterations. Default = ``1000``"""
    guess: Guess = field(default_factory=lambda: Guess())
    """Initial guess. Default = ``stls.Guess()``"""
    theory: str = "STLS2D"


@serialize.serializable_dataclass
class Result(hf2d.Result):
    """
    Class used to store the results for the :obj:`qupled.stls2d.Stls2D` class.
    """

    error: float = None
    """Final error of the scheme. Default = ``None``"""


@serialize.serializable_dataclass
class Guess:
    wvg: np.ndarray = None
    """Wave-vector grid. Default = ``None``"""
    ssf: np.ndarray = None
    """Static structure factor. Default = ``None``"""

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


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
