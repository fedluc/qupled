#!/usr/bin/env python

import sys
import qupled.qupled as qp
import qupled.Input as Input

class Stls():
            
    # Constructor
    def __init__(self,
                 coupling,
                 degeneracy,
                 chemicalPotential = None,
                 cutoff = None,
                 error = None,
                 mixing = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 resolution = None ):
        # Allowed theories
        self.allowedTheories = ["STLS"]
        # Default inputs
        self.inputs = Input.Stls(coupling, degeneracy, "STLS")
        # Scheme to solve and associated input
        self.scheme = None
        self.__schemeInputs = None
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
        
    # Compute
    def compute(self):
        self.checkInputs()
        self.__schemeInputs = self.inputs
        self.scheme = qp.Stls(self.__schemeInputs.inputs)
        self.scheme.compute()

    # Check input before computing
    def checkInputs(self):
        if (self.inputs.getTheory() not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")
        
class StlsIet(Stls):
            
    # Constructor
    def __init__(self,
                 coupling,
                 degeneracy,
                 theory,
                 chemicalPotential = None,
                 cutoff = None,
                 error = None,
                 mapping = None,
                 mixing = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 scheme2DIntegrals = None,
                 resolution = None ):
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, iterations, matsubara,
                         outputFrequency, restartFile,
                         resolution)
        # Allowed theories
        self.allowedTheories = ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        # Set theory
        self.inputs.setTheory(theory)
        self.checkInputs()
        # Non-default inputs
        if (mapping != None):
            self.inputs.IETMapping(iet)
        if (scheme2DIntegrals != None):
            self.inputs.setInt2DScheme(scheme2DIntegrals)
