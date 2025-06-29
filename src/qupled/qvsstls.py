from __future__ import annotations

from . import native
from . import qstls
from . import vsstls
from . import serialize


class Solver(vsstls.Solver):
    """
    Class used to solve the QVStls scheme.
    """

    # Native classes used to solve the scheme
    native_scheme_cls = native.QVSStls
    native_inputs_cls = native.QVSStlsInput

    def __init__(self):
        super().__init__()
        self.results: vsstls.Result = vsstls.Result()

    def compute(self, inputs: Input):
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        qstls.Solver.find_fixed_adr_in_database(self, inputs)
        super().compute(inputs)

    def _update_input_data(self, inputs: Input):
        """
        Updates the input data with additional attributes specific to the current instance.

        This method overrides the parent class's `_update_input_data` method to include
        logic for setting a default `fixed_run_id` if it is not already provided in the
        `inputs`.

        Args:
            inputs: Input parameters.

        Side Effects:
            - If `inputs.fixed_run_id` is `None`, it is set to the value of `self.run_id`.
        """
        super()._update_input_data(inputs)
        if inputs.fixed_run_id is None:
            inputs.fixed_run_id = self.run_id


@serialize.serializable_dataclass
class Input(vsstls.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qvsstls.QVSStls` class.
    """

    theory: str = "QVSSTLS"


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, vsstls.Result)
