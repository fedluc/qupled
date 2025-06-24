from __future__ import annotations

import numpy as np

from . import native
from . import output
from . import stls


class StlsIet(stls.Stls):
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


class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.stlsiet.StlsIet` class.
    Accepted theories: ``STLS-HNC``, ``STLS-IOI`` and ``STLS-LCT``.
    """

    def __init__(self, coupling: float, degeneracy: float, theory: str):
        super().__init__(coupling, degeneracy)
        if theory not in {"STLS-HNC", "STLS-IOI", "STLS-LCT"}:
            raise ValueError("Invalid dielectric theory")
        self.theory = theory
        self.mapping = "standard"
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
        self.guess: Guess = Guess()
        """Initial guess. Default = ``stlsiet.Guess()``"""


class Result(stls.Result):
    """
    Class used to store the results for the :obj:`qupled.stlsiet.StlsIet` class.
    """

    def __init__(self):
        super().__init__()
        self.bf: np.ndarray = None
        """Bridge function adder"""


class Guess(stls.Guess):

    def __init__(
        self,
        wvg: np.ndarray = None,
        ssf: np.ndarray = None,
        lfc: np.ndarray = None,
    ):
        super().__init__(wvg, ssf)
        self.lfc = lfc
        """ Local field correction. Default = ``None``"""


if __name__ == "__main__":
    from .mpi_worker import run_mpi_worker

    run_mpi_worker(Input, Result, StlsIet.native_inputs_cls, StlsIet.native_scheme_cls)
