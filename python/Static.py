#!/usr/bin/env python

import sys
import os
import tarfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import qupled.qupled as qp

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
        if (chemicalPotential is not None): self.inputs.chemicalPotential = chemicalPotential
        if (cutoff is not None): self.inputs.cutoff = cutoff
        if (error is not None): self.inputs.error = error
        if (mixing is not None): self.inputs.mixing = mixing
        if (iterations is not None): self.inputs.iterations = iterations
        if (matsubara is not None): self.inputs.matsubara = matsubara
        if (outputFrequency is not None): self.inputs.outputFrequency = outputFrequency
        if (restartFile is not None): self.inputs.restartFile = restartFile
        if (resolution is not None): self.inputs.resolution = resolution

    # Compute
    def compute(self):
        self.checkInputs()
        self.__schemeInputs = self.inputs
        self.scheme = qp.Stls(self.__schemeInputs)
        self.scheme.compute()
        self.save()
        
    # Check input before computing
    def checkInputs(self):
        if (self.inputs.theory not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")
      
    # Plot results        
    def plot(self, toPlot, matsubaraFrequencies = None, rdfGrid = None):
        wvg = self.scheme.wvg
        xlabel = "x = k/$k_F$"
        if ("idr" in toPlot): self.plotIdr(matsubaraFrequencies)
        if ("rdf" in toPlot): self.plotRdf(rdfGrid)
        if ("sdr" in toPlot): self.plot1D(wvg, self.scheme.sdr, xlabel, "$\chi$(x)")
        if ("slfc" in toPlot): self.plot1D(wvg, self.scheme.slfc, xlabel, "G(x)")
        if ("ssf" in toPlot): self.plot1D(wvg, self.scheme.ssf, xlabel, "S(x)")
        if ("ssfHF" in toPlot): self.plot1D(wvg, self.scheme.ssfHF, xlabel, "$S_{HF}$(x)")
        
    def plotIdr(self, matsubaraFrequencies):
        if (matsubaraFrequencies is None) : matsubaraFrequencies = np.arange(self.inputs.matsubara)
        self.plot1DParametric(self.scheme.wvg, self.scheme.idr,
                              "x = k/k_F", "$\phi$(x)",
                              matsubaraFrequencies)
        
    def plotRdf(self, rdfGrid):
        if (rdfGrid is None) : rdfGrid = np.arange(0.01, 10.0, 0.01)
        rdf = self.scheme.getRdf(rdfGrid)
        self.plot1D(rdfGrid, rdf, "x = rk_F", "g(r)")

    def plot1D(self, x, y, xlabel, ylabel):
        plt.plot(x, y, "b")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()        

    def plot1DParametric(self, x, y, xlabel, ylabel, parameters):
        numParameters = parameters.size
        cmap = cm.get_cmap(name="viridis")
        for i in np.arange(numParameters):
            color = cmap(1.0*i/numParameters)
            plt.plot(x, y[:,parameters[i]], color=color)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()        

    def save(self):
        inputs = pd.DataFrame({
            "coupling" : self.__schemeInputs.coupling,
            "degeneracy" : self.__schemeInputs.degeneracy,
            "theory" : self.__schemeInputs.theory,
            "error" : self.__schemeInputs.error,
            "resolution" : self.__schemeInputs.resolution,
            "cutoff" : self.__schemeInputs.cutoff,
            "matsubara" : self.__schemeInputs.matsubara
            }, index=["inputs"]).to_csv("inputs.csv")
        pd.DataFrame(self.scheme.idr).to_csv("idr.csv")
        pd.DataFrame(self.scheme.sdr).to_csv("sdr.csv")
        pd.DataFrame(self.scheme.slfc).to_csv("slfc.csv")
        pd.DataFrame(self.scheme.ssf).to_csv("ssf.csv")
        pd.DataFrame(self.scheme.ssfHF).to_csv("ssfHF.csv")
        pd.DataFrame({"uint" : self.scheme.uInt}, index=[0]).to_csv("uint.csv")
        pd.DataFrame(self.scheme.wvg).to_csv("wvg.csv")
        file = tarfile.open("data.tar.gz","w:gz")
        file.add("inputs.csv")
        file.add("idr.csv")
        file.add("sdr.csv")
        file.add("slfc.csv")
        file.add("ssf.csv")
        file.add("ssfHF.csv")
        file.add("uint.csv")
        file.add("wvg.csv")
        file.close()
        os.remove("inputs.csv")
        os.remove("idr.csv")
        os.remove("sdr.csv")
        os.remove("slfc.csv")
        os.remove("ssf.csv")
        os.remove("ssfHF.csv")
        os.remove("uint.csv")
        os.remove("wvg.csv")
        
        
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
        if (mapping is not None): self.inputs.iet = mapping
        if (scheme2DIntegrals is not None): self.inputs.int2DScheme = scheme2DIntegrals
            
    # Plot results        
    def plot(self, toPlot, extraParameters = None):
        super().plot(toPlot, extraParameters)
        if ("bf" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.bf, "x = k/k_F", "B(x)")
