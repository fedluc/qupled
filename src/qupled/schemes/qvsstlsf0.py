from __future__ import annotations

from qupled import native
from qupled.database.scheme_tables import BaseTableKeys, TableKeys
from qupled.schemes import qvsstls, vsstls
from qupled.util import serialize


class Solver(vsstls.Solver):
    """
    Class used to solve the QVSStls-F0 scheme.
    """

    native_scheme_cls = native.QVSStlsF0
    native_inputs_cls = native.QVSStlsF0Input

    def __init__(self):
        super().__init__()
        self.results: vsstls.Result = vsstls.Result()

    @staticmethod
    def _same_optional_float(a, b, atol: float = 1.0e-14) -> bool:
        if a is None and b is None:
            return True
        if a is None or b is None:
            return False
        return abs(float(a) - float(b)) <= atol

    def compute(self, inputs: Input):
        self.find_fixed_adr_in_database(inputs)
        super().compute(inputs)

    def find_fixed_adr_in_database(self, inputs: "Input"):
        scheme_tables = self.db_handler.scheme_tables
        runs = scheme_tables.inspect_runs()
        inputs.fixed_run_id = None
        for run in runs:
            same_degeneracy = run[TableKeys.DEGENERACY.value] == inputs.degeneracy
            same_theory = run[TableKeys.THEORY.value] == inputs.theory
            same_coupling = run[TableKeys.COUPLING.value] == inputs.coupling
            if not (same_theory and same_degeneracy and same_coupling):
                continue
            run_id = run[BaseTableKeys.PRIMARY_KEY.value]
            run_inputs = scheme_tables.get_inputs(run_id)
            if (
                run_inputs["cutoff"] == inputs.cutoff
                and run_inputs["matsubara"] == inputs.matsubara
                and run_inputs["resolution"] == inputs.resolution
                and self._same_optional_float(
                    run_inputs.get("pimc_eta"), inputs.pimc_eta
                )
                and self._same_optional_float(
                    run_inputs.get("pimc_y_sec"), inputs.pimc_y_sec
                )
                and self._same_optional_float(
                    run_inputs.get("pimc_a_cutoff"), inputs.pimc_a_cutoff
                )
            ):
                print(f"Loading fixed ADR from database for run_id = {run_id}")
                inputs.fixed_run_id = run_id
                return

    def _update_input_data(self, inputs: Input):
        super()._update_input_data(inputs)
        if inputs.fixed_run_id is None:
            inputs.fixed_run_id = self.run_id


@serialize.serializable_dataclass
class Input(qvsstls.Input):
    """
    Class used to manage the input for the :obj:`qupled.qvsstlsf0.QVSStlsF0` class.
    """

    theory: str = "QVSSTLS-F0"
    pimc_eta: float | None = None
    pimc_y_sec: float | None = None
    pimc_a_cutoff: float | None = None


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, vsstls.Result)
