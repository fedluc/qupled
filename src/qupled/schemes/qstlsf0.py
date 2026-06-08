from __future__ import annotations

from qupled import native
from qupled.database.scheme_tables import BaseTableKeys, TableKeys
from qupled.schemes import qstls, stls
from qupled.util import serialize


class Solver(qstls.Solver):
    """
    Class used to solve the QSTLS-F0 scheme (FD + asymptotic tail).
    """

    native_scheme_cls = native.QstlsF0
    native_inputs_cls = native.QstlsPimcInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()

    @staticmethod
    def _same_optional_float(a, b, atol: float = 1.0e-14) -> bool:
        if a is None and b is None:
            return True
        if a is None or b is None:
            return False
        return abs(float(a) - float(b)) <= atol

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


@serialize.serializable_dataclass
class Input(qstls.Input):
    """
    Class used to manage the input for :obj:`qupled.qstlsf0.QstlsF0`.
    """

    theory: str = "QSTLS-F0"
    pimc_eta: float | None = None
    pimc_y_sec: float | None = None
    pimc_a_cutoff: float | None = None


@serialize.serializable_dataclass
class Result(stls.Result):
    """
    Results for :obj:`qupled.qstlsf0.QstlsF0`.
    """

    f0_grid: any = None
    f0_values: any = None


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
