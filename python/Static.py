#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
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
        self.inputs = qp.StlsInput(coupling, degeneracy, "STLS")
        # Scheme to solve and associated input and solution
        self.scheme = None
        self.__schemeInputs = None
        self.__schemeSolution = None     
        # Non-default inputs
        if (chemicalPotential != None):
            sys.exit("Setting of chemical potential guess is temporary disabled")
            # self.inputs.setChemicalPotential(-10, 10)
        if (cutoff != None): self.inputs.cutoff = cutoff
        if (error != None): self.inputs.error = error
        if (mixing != None): self.inputs.mixing = mixing
        if (iterations != None): self.inputs.iterations = iterations
        if (matsubara != None): self.inputs.matsubara = matsubara
        if (outputFrequency != None): self.inputs.outputFrequency = outputFrequency
        if (restartFile != None): self.inputs.restartFile = restartFile
        if (resolution != None): self.inputs.resolution = resolution

    # Compute
    def compute(self):
        self.checkInputs()
        self.__schemeInputs = self.inputs
        self.scheme = qp.Stls(self.__schemeInputs)
        self.scheme.compute()
        # self.getSolution() # Otherwise memory is cleaned after c++ is done
        
    # Check input before computing
    def checkInputs(self):
        if (self.inputs.theory not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")
      
    # Plot results        
    def plot(self, toPlot, extraParameters = None):
        if ("idr" in toPlot): sys.exit("Plot of ideal density response is not yet available")
        if ("rdf" in toPlot): sys.exit("Plot of radial distribution function is not yet available")
        if ("sdr" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.sdr, "b", "x = k/k_F", "$\chi$(x)")
        if ("slfc" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.slfc, "b", "x = k/k_F", "G(x)")
        if ("ssf" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.ssf, "b", "x = k/k_F", "S(x)")
        if ("ssfHF" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.ssfHF, "b", "x = k/k_F", "$S_{HF}$(x)")

    def plot1D(self, x, y, color, xlabel, ylabel):
        plt.plot(x, y, color)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()        

             
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
        self.inputs.theory = theory
        self.checkInputs()
        # Non-default inputs
        if (mapping != None): self.inputs.iet = mapping
        if (scheme2DIntegrals != None): self.inputs.int2DScheme = scheme2DIntegrals
            
    # Plot results        
    def plot(self, toPlot, extraParameters = None):
        super().plot(toPlot, extraParameters)
        if ("bf" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.bf, "b", "x = k/k_F", "B(x)")
