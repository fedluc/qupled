from __future__ import annotations

from dataclasses import field
import numpy as np

from . import native
from . import output
from . import stls
from . import serialize


class Solver(stls.Solver):
    """
    Class used to solve the StlsIet schemes.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.StlsIet
    native_inputs_cls = native.StlsIetInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()

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
        names = ["wvg", "ssf", "lfc"]
        results = output.DataBase.read_results(run_id, database_name, names)
        return Guess(
            results[names[0]],
            results[names[1]],
            results[names[2]],
        )


@serialize.serializable_dataclass
class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.stlsiet.StlsIet` class.
    Accepted theories: ``STLS-HNC``, ``STLS-IOI`` and ``STLS-LCT``.
    """

    mapping: str = "standard"
    r"""
        Mapping for the classical-to-quantum coupling parameter
        :math:`\Gamma` used in the iet schemes. Allowed options include:

        - standard: :math:`\Gamma \propto \Theta^{-1}`

        - sqrt: :math:`\Gamma \propto (1 + \Theta)^{-1/2}`

        - linear: :math:`\Gamma \propto (1 + \Theta)^{-1}`

        where :math:`\Theta` is the degeneracy parameter. Far from the ground state
        (i.e. :math:`\Theta\gg1`) all mappings lead identical results, but at
        the ground state they can differ significantly (the standard
        mapping diverges). Default = ``standard``.
        """
    guess: Guess = field(default_factory=lambda: Guess())

    def __post_init__(self):
        if self.theory not in {"STLS-HNC", "STLS-IOI", "STLS-LCT"}:
            raise ValueError("Invalid dielectric theory")


@serialize.serializable_dataclass
class Result(stls.Result):
    """
    Class used to store the results for the :obj:`qupled.stlsiet.StlsIet` class.
    """

    bf: np.ndarray = None
    """Bridge function adder"""


@serialize.serializable_dataclass
class Guess(stls.Guess):
    lfc: np.ndarray = None
    """ Local field correction. Default = ``None``"""


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, Result, Solver.native_inputs_cls, Solver.native_scheme_cls)
