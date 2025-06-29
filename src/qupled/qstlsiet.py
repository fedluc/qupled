from __future__ import annotations

from . import native
from . import qstls
from . import serialize
from . import stlsiet


class Solver(qstls.Solver):
    """
    Class used to solve the Qstls-IET schemes.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.QstlsIet
    native_inputs_cls = native.QstlsIetInput

    def __init__(self):
        super().__init__()
        self.results: stlsiet.Result = stlsiet.Result()

    @staticmethod
    def get_initial_guess(
        run_id: str, database_name: str | None = None
    ) -> stlsiet.Guess:
        return stlsiet.Solver.get_initial_guess(run_id, database_name)


@serialize.serializable_dataclass
class Input(stlsiet.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qstlsiet.QStlsIet` class.
    Accepted theories: ``QSTLS-HNC``, ``QSTLS-IOI`` and ``QSTLS-LCT``.
    """

    allowed_theories = {"QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"}


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, stlsiet.Result)
