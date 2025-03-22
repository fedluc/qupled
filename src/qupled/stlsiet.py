# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------

from __future__ import annotations
import sys
import pandas as pd
import numpy as np
from . import native
from . import util
from . import rpa
from . import stls
from . import base


class StlsIet(base.IterativeScheme):

    # Compute
    @util.MPI.record_time
    @util.MPI.synchronize_ranks
    def compute(self, inputs: StlsIet.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        scheme = native.Stls(inputs.to_native())
        self._compute(scheme)
        self._save(scheme, inputs)

    # Save results to disk
    @util.MPI.run_only_on_root
    def _save(self, scheme, inputs) -> None:
        """Stores the results obtained by solving the scheme."""
        super()._save(scheme, inputs, self.Results(scheme))
        pd.DataFrame(scheme.bf).to_hdf(self.hdf_file_name, key="bf")

    # Input class
    class Input(stls.Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.StlsIet` class.
        Accepted theories: ``STLS-HNC``, ``STLS-IOI`` and ``STLS-LCT``.
        """

        def __init__(self, coupling: float, degeneracy: float, theory: str):
            super().__init__(coupling, degeneracy)
            if theory not in {"STLS-HNC", "STLS-IOI", "STLS-LCT"}:
                sys.exit("Invalid dielectric theory")
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

    # Results class
    class Results(rpa.Rpa.Results):
        """
        Class used to store the results for the :obj:`qupled.classic.StlsIet` class.
        """

        def __init__(self, scheme):
            super().__init__(scheme, init_from_native=False)
            self.bf: np.ndarray = None
            """Bridge function adder"""
            base.Result.from_native(self, scheme)
