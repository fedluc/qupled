from __future__ import annotations

from . import native
from . import qstls
from . import vsstls


class QVSStls(vsstls.VSStls):
    """
    Class used to solve the QVStls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: Result = Result()
        # Undocumented properties
        self.native_scheme_cls = native.QVSStls
        self.native_inputs = native.QVSStlsInput()

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

    @staticmethod
    def get_initial_guess(run_id: int, database_name: str | None = None) -> qstls.Guess:
        return qstls.Qstls.get_initial_guess(run_id, database_name)


# Input class
class Input(vsstls.Input, qstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qvsstls.QVSStls` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        vsstls.Input.__init__(self, coupling, degeneracy)
        qstls.Input.__init__(self, coupling, degeneracy)
        # Undocumented default values
        self.theory: str = "QVSSTLS"


# Result
class Result(vsstls.Result, qstls.Result):
    """
    Class used to manage the results for the :obj:`qupled.qvsstls.QVSStls` class.
    """

    def __init__(self):
        super().__init__()
