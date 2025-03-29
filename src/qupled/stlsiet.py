# -----------------------------------------------------------------------
# StlsIet class
# -----------------------------------------------------------------------

from __future__ import annotations
import sys
import pandas as pd
import numpy as np
from . import native
from . import stls
from . import base


class StlsIet(base.IterativeScheme):

    # Compute
    def compute(self, inputs: StlsIet.Input) -> None:
        """
        Solves the scheme and saves the results.

        Args:
            inputs: Input parameters.
        """
        super().compute(inputs, native.Stls, native.StlsInput(), self.Result())

    # Input class
    class Input(stls.Stls.Input):
        """
        Class used to manage the input for the :obj:`qupled.classic.StlsIet` class.
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

    # Results class
    class Result(stls.Stls.Result):
        """
        Class used to store the results for the :obj:`qupled.classic.StlsIet` class.
        """

        def __init__(self):
            super().__init__()
            self.bf: np.ndarray = None
            """Bridge function adder"""
