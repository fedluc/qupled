#!/usr/bin/env python

import sys
import qupled.qupled as qp
import qupled.Input as Input

class Stls():

    # Constructor
    def __init__(self,
                 coupling,
                 degeneracy
                 chemicalPotential = None,
                 cutoff = None,
                 error = None,
                 mixing = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 resolution = None ):
        # Default inputs
        self.inputs = Input.StlsInput(coupling, degeneracy, "STLS")
        # Non-default inputs
        if (chemicalPotential != None):
            self.inputs.setChemicalPotential(-10, 10)
        if (cutoff != None):
            self.inputs.setWaveVectorGridCutoff(cutoff)
        if (error != None):
            self.inputs.setErrMin(error)
        if (mixing != None):
            self.inputs.setMixingParameter(mixing)
        if (iterations != None):
            self.inputs.setNIter(iterations)
        if (matsubara != None):
            self.inputs.setNMatsubara(matsubara)
        if (outputFrequency != None):
            self.inputs.setOutIter(outputFrequency)
        if (restartFile != None):
            self.inputs.setRestartFile(restartFile)
        if (resolution != None):
            self.inputs.setWaveVectorGridRes(resolution)
        # Scheme to solve and associated input
        self.scheme = None
        self.__schemeInputs = None
        
    # Compute
    def compute(self):
        self.__schemeInputs = self.inputs
        self.scheme = qp.Stls(self.__schemeInputs.inputs)
        self.scheme.compute()


class StlsIet(Stls):

    # Constructor
    def __init__(self,
                 coupling,
                 degeneracy,
                 theory,
                 chemicalPotential = None,
                 cutoff = None,
                 error = None,
                 iet = None,
                 mixing = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 scheme2DIntegrals = None
                 resolution = None ):
        # Call parent constructor
        super().__init___(self, coupling, degeneracy,
                          chemicalPotential, cutoff, error,
                          mixing, iterations, matsubara,
                          outputFrequency, restartFile,
                          resolution)
        # Set theory
        if (theory == "STLS-HNC" or
            theory == "STLS-IOI" or
            theory == "STLS-LCT"):
            self.inputs.setTheort(theory)
        else:
            sys.exit("StlsIet: Unknown theory")
        # Non-default inputs
        if (iet != None):
            self.inputs.IETMapping(iet)
        if (scheme2DIntegrals != None):
            self.inputs.setInt2DScheme(scheme2DIntegrals)
