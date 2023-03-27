#!/usr/bin/env python

import sys
import qupled.qupled as qp
import qupled.Input as Input

class Stls():

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
                 scheme2DIntegrals = None,
                 threads = None,
                 resolution = None ):
        # Default inputs
        self.inputs = Input.StlsInput(coupling, degeneracy, theory)
        # Non-default inputs
        if (chemicalPotential != None):
            self.inputs.setChemicalPotential(-10, 10)
        if (cutoff != None):
            self.inputs.setWaveVectorGridCutoff(cutoff)
        if (error != None):
            self.inputs.setErrMin(error)
        if (iet != None):
            self.inputs.IETMapping(iet)
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
        if (scheme2DIntegrals != None):
            self.inputs.setInt2DScheme(scheme2DIntegrals)
        if (threads != None):
            self.inputs.setNThreads(threads)
        if (resolution != None):
            self.inputs.setWaveVectorGridRes(resolution)
        
    # Compute
    def compute(self):
        # Create something like self.inputs in latest compute
        self.scheme = qp.Stls(self.inputs.inputs)
        self.scheme.compute()
