#!/usr/bin/env python

import sys
import os
import glob
import zipfile as zf
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
                 guess = None,
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
        self.schemeInputs = None
        self.schemeSolution = None
        # File to store output on disk
        self.hdfFileName = None
        # Non-default inputs
        if (chemicalPotential is not None): self.inputs.chemicalPotential = chemicalPotential
        if (cutoff is not None): self.inputs.cutoff = cutoff
        if (error is not None): self.inputs.error = error
        if (mixing is not None): self.inputs.mixing = mixing
        if (guess is not None): self.inputs.guess = guess
        if (iterations is not None): self.inputs.iterations = iterations
        if (matsubara is not None): self.inputs.matsubara = matsubara
        if (outputFrequency is not None): self.inputs.outputFrequency = outputFrequency
        if (restartFile is not None): self.inputs.restartFile = restartFile
        if (resolution is not None): self.inputs.resolution = resolution

    # Compute
    def compute(self):
        self.checkInputs()
        self.schemeInputs = self.inputs
        self.scheme = qp.Stls(self.schemeInputs)
        self.scheme.compute()
        self.setHdfFile()
        self.save()
        
    # Check input before computing
    def checkInputs(self):
        if (self.inputs.theory not in self.allowedTheories):
            sys.exit("Invalid dielectric theory")

    # Compute radial distribution function
    def computeRdf(self, rdfGrid, writeToHdf = True):
        if (self.schemeInputs == None):
            sys.exit("No solution to compute the radial distribution function")
        rdf = self.scheme.getRdf(rdfGrid)
        if (writeToHdf) : pd.DataFrame(rdf).to_hdf(self.hdfFileName, key="rdf", mode="r+")
        return rdf
        
    # Plot results        
    def plot(self, toPlot, matsubara = None, rdfGrid = None):
        wvg = self.scheme.wvg
        xlabel = "Wave vector"
        if ("idr" in toPlot): self.plotIdr(matsubara)
        if ("rdf" in toPlot): self.plotRdf(rdfGrid)
        if ("sdr" in toPlot): self.plot1D(wvg, self.scheme.sdr, xlabel, "Static density response")
        if ("slfc" in toPlot): self.plot1D(wvg, self.scheme.slfc, xlabel, "Static local field correction")
        if ("ssf" in toPlot): self.plot1D(wvg, self.scheme.ssf, xlabel, "Static structure factor")
        if ("ssfHF" in toPlot): self.plot1D(wvg, self.scheme.ssfHF, xlabel, "Hartree-Fock static structure factor")
        
    def plotIdr(self, matsubara):
        if (matsubara is None) : matsubara = np.arange(self.inputs.matsubara)
        self.plot1DParametric(self.scheme.wvg, self.scheme.idr,
                              "Wave vector", "Ideal density response",
                              matsubara)
        
    def plotRdf(self, rdfGrid):
        if (rdfGrid is None) : rdfGrid = np.arange(0.01, 10.0, 0.01)
        self.plot1D(rdfGrid, self.computeRdf(rdfGrid), "Inter-particle distance", "radial distribution function")

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

    # Save results to disk
    def setHdfFile(self):
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.schemeInputs.coupling,
                                                           self.schemeInputs.degeneracy,
                                                           self.schemeInputs.theory)
    
    def save(self):
        pd.DataFrame({
            "coupling" : self.schemeInputs.coupling,
            "degeneracy" : self.schemeInputs.degeneracy,
            "theory" : self.schemeInputs.theory,
            "error" : self.schemeInputs.error,
            "resolution" : self.schemeInputs.resolution,
            "cutoff" : self.schemeInputs.cutoff,
            "matsubara" : self.schemeInputs.matsubara
            }, index=["inputs"]).to_hdf(self.hdfFileName, key="inputs", mode="w")
        pd.DataFrame(self.scheme.idr).to_hdf(self.hdfFileName, key="idr")
        pd.DataFrame(self.scheme.sdr).to_hdf(self.hdfFileName, key="sdr")
        pd.DataFrame(self.scheme.slfc).to_hdf(self.hdfFileName, key="slfc")
        pd.DataFrame(self.scheme.ssf).to_hdf(self.hdfFileName, key="ssf")
        pd.DataFrame(self.scheme.ssfHF).to_hdf(self.hdfFileName, key="ssfHF")
        pd.DataFrame({"uint" : self.scheme.uInt}, index=[0]).to_hdf(self.hdfFileName, key="uint")
        pd.DataFrame(self.scheme.wvg).to_hdf(self.hdfFileName, key="wvg")
        
        
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
                 guess = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 scheme2DIntegrals = None,
                 resolution = None ):
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, guess, iterations,
                         matsubara, outputFrequency,
                         restartFile, resolution)
        # Allowed theories
        self.allowedTheories = ["STLS-HNC", "STLS-IOI", "STLS-LCT"]
        # Set theory
        self.inputs.theory = theory
        self.checkInputs()
        # File to store output on disk
        self.__hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                           self.inputs.degeneracy,
                                                           self.inputs.theory)
        # Non-default inputs
        if (mapping is not None): self.inputs.iet = mapping
        if (scheme2DIntegrals is not None): self.inputs.int2DScheme = scheme2DIntegrals
            
    # Plot results        
    def plot(self, toPlot, matsubara = None, rdfGrid = None):
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        
    # Save results to disk
    def save(self):
        super().save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.__hdfFileName, key="bf")



class Qstls(Stls):
            
    # Constructor
    def __init__(self,
                 coupling,
                 degeneracy,
                 chemicalPotential = None,
                 cutoff = None,
                 error = None,
                 mixing = None,
                 guess = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 resolution = None,
                 threads = None):
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, guess, iterations,
                         matsubara, outputFrequency,
                         restartFile, resolution)
        # Allowed theories
        self.allowedTheories = ["QSTLS"]
        # Set theory
        self.inputs.theory = "QSTLS"
        # Qstls inputs
        self.qInputs = qp.QstlsInput()
        self.schemeqInputs = None;
        self.checkInputs()
        # File to store output on disk
        self.hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                         self.inputs.degeneracy,
                                                         self.inputs.theory)
        # Non-default inputs
        if (threads is not None): self.inputs.threads = threads

    # Compute
    def compute(self):
        self.checkInputs()
        self.schemeInputs = self.inputs
        self.schemeqInputs = self.qInputs
        self.scheme = qp.Qstls(self.schemeInputs, self.schemeqInputs)
        self.scheme.compute()
        self.setHdfFile()
        self.save()
        
    # Plot results        
    def plot(self, toPlot, matsubara = None, rdfGrid = None):
        super().plot(toPlot, matsubara, rdfGrid)
        if ("adr" in toPlot): self.plotAdr(matsubara)

    def plotAdr(self, matsubara):
        if (matsubara is None) : matsubara = np.arange(self.inputs.matsubara)
        self.plot1DParametric(self.scheme.wvg, self.scheme.adr,
                              "Wave vector", "Auxiliary density response",
                              matsubara)
        
    # Save results to disk
    def save(self):
        super().save()
        pd.DataFrame(self.scheme.adr).to_hdf(self.hdfFileName, key="adr")
        shape = self.scheme.adr_fixed.shape
        pd.DataFrame(self.scheme.adr_fixed.reshape(shape[0], shape[1]*shape[2])).to_hdf(self.hdfFileName, key="adr_fixed")

        

class QstlsIet(Qstls):
            
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
                 guess = None,
                 iterations = None,
                 matsubara = None,
                 outputFrequency = None,
                 restartFile = None,
                 scheme2DIntegrals = None,
                 resolution = None,
                 threads = None):
        # Call parent constructor
        super().__init__(coupling, degeneracy,
                         chemicalPotential, cutoff, error,
                         mixing, guess, iterations,
                         matsubara, outputFrequency,
                         restartFile, resolution)
        # Allowed theories
        self.allowedTheories = ["QSTLS-HNC", "QSTLS-IOI", "QSTLS-LCT"]
        # Set theory
        self.inputs.theory = theory
        self.checkInputs()
        # File to store output on disk
        self.__hdfFileName = "rs%5.3f_theta%5.3f_%s.h5" % (self.inputs.coupling,
                                                           self.inputs.degeneracy,
                                                           self.inputs.theory)
        # Non-default inputs
        if (mapping is not None): self.inputs.iet = mapping
        if (scheme2DIntegrals is not None): self.inputs.int2DScheme = scheme2DIntegrals
            
    # Plot results        
    def plot(self, toPlot, matsubara = None, rdfGrid = None):
        super().plot(toPlot, matsubara, rdfGrid)
        if ("bf" in toPlot): self.plot1D(self.scheme.wvg, self.scheme.bf, "Wave vector", "Bridge function adder")
        
    # Save results to disk
    def save(self):
        super().save()
        pd.DataFrame(self.scheme.bf).to_hdf(self.__hdfFileName, key="bf")
        # Zip all files for the fixed component of the auxiliary density response
        adrFileName = "adr_iet_fixed_rs%5.3f_theta%5.3f_%s" % (self.inputs.coupling,
                                                               self.inputs.degeneracy,
                                                               self.inputs.theory)
        zipFile = zf.ZipFile(adrFileName + ".zip", "w")
        for adrFile in glob.glob(adrFileName + "_wv*.bin"):
            zipFile.write(adrFile)
            os.remove(adrFile)
