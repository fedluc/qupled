from __future__ import annotations

from . import database
from . import native
from . import stls


class Qstls(stls.Stls):
    """
    Class used to solve the Qstls scheme.
    """

    def __init__(self):
        super().__init__()
        self.results: stls.Result = stls.Result()
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


# Input class
class Input(stls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qstls.Qstls` class.
    """

    def __init__(self, coupling: float, degeneracy: float):
        super().__init__(coupling, degeneracy)
        # Undocumented default values
        self.fixed_run_id: int | None = None
        self.theory = "QSTLS"
