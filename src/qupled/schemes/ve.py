from __future__ import annotations

from dataclasses import field
import numpy as np

from qupled import native
from qupled.schemes import qstls, stls, vsstls
from qupled.util import serialize


class Solver(qstls.Solver):
    """
    Class used to solve the viscoelastic qSTLS scheme.
    """

    native_scheme_cls = native.VE
    native_inputs_cls = native.VEInput

    def __init__(self):
        super().__init__()
        self.results: Result = Result()

    def compute(self, inputs: "Input"):
        qstls.Solver.find_fixed_adr_in_database(self, inputs)
        self._fill_free_energy_integrand(inputs)
        super().compute(inputs)

    def _fill_free_energy_integrand(self, inputs: "Input"):
        target_coupling = inputs.coupling
        missing_state_points = self._get_missing_state_points(inputs)
        did_subcalls = len(missing_state_points) > 0
        for coupling in missing_state_points:
            print("---------------------------------------------------------------")
            print(f"Subcall: solving {inputs.theory} scheme for rs = {coupling:.5f}")
            inputs.coupling = coupling
            self.compute(inputs)
            self._update_input_data(inputs)
        if did_subcalls:
            print("---------------------------------------------------------------")
            print("Subcalls completed.")
        inputs.coupling = target_coupling

    @staticmethod
    def _get_missing_state_points(inputs: "Input") -> np.ndarray:
        rs = inputs.coupling
        drs = inputs.coupling_resolution
        expected_grid = np.arange(drs, rs - 0.1 * drs, drs)
        actual_grid = inputs.free_energy_integrand.grid
        precision = int(np.abs(np.log10(drs))) if drs < 1.0 else 0
        return (
            np.setdiff1d(
                np.round(expected_grid, precision), np.round(actual_grid, precision)
            )
            if actual_grid is not None
            else expected_grid
        )

    def _update_input_data(self, inputs: "Input"):
        free_energy_integrand = vsstls.FreeEnergyIntegrand(
            self.results.free_energy_grid, self.results.free_energy_integrand
        )
        inputs.free_energy_integrand = free_energy_integrand
        if inputs.fixed_run_id is None:
            inputs.fixed_run_id = self.run_id


@serialize.serializable_dataclass
class Input(qstls.Input):
    """
    Class used to manage the input for :obj:`qupled.ve.VE`.
    """

    theory: str = "VE"
    coupling_resolution: float = 0.1
    free_energy_integrand: vsstls.FreeEnergyIntegrand = field(
        default_factory=lambda: vsstls.FreeEnergyIntegrand()
    )


@serialize.serializable_dataclass
class Result(stls.Result):
    """
    Results for :obj:`qupled.ve.VE`.
    """

    free_energy_grid: np.ndarray = None
    free_energy_integrand: np.ndarray = None
    matsubara_grid: np.ndarray = None
    a_coeff: np.ndarray = None
    a1_coeff: np.ndarray = None
    mxc_l: np.ndarray = None
    a0_coeff: float = None
    kxc0: float = None
    mxc_inf: float = None
    cxc: float = None
    omega_m2: float = None


if __name__ == "__main__":
    Solver.run_mpi_worker(Input, Result)
