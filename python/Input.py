#!/usr/bin/env python

import sys
import qupled.qupled as qp

class Input:

    # Constructor
    def __init__(self, coupling, degeneracy, theory):
        goodInput = ( isinstance(coupling, (int, float)) and \
                      isinstance(degeneracy, (int, float)) and \
                      isinstance(theory, str) )
        if (goodInput):
            self.inputs = qp.Input(coupling, degeneracy, theory)
        else:
            self.throwInputFormatError("Input:")
        self.isDefault = True
            
    # Set number of OMP threads
    def setNThreads(self, num):
        if (isinstance(num, int)):
            self.inputs.setNThreads(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setNThreads")
    
    # Set scheme for 2D integrals
    def setInt2DScheme(self, name):
        if (isinstance(num, str)):
            self.inputs.setInt2DScheme(name)
            self.isDefault = False
        else:
            self.throwInputFormatError("setInt2DScheme")

    # Print content of input structure
    def print(self):
        self.inputs.print()
            

    # Throw error messages
    def throwInputFormatError(self, callerName):
        sys.exit(callerName + ": Wrong input format")


class StlsInput(Input):

    # Constructor
    def __init__(self, coupling, degeneracy, theory):
        goodInput = ( isinstance(coupling, (int, float)) and \
                      isinstance(degeneracy, (int, float)) and \
                      isinstance(theory, str) )
        if (goodInput):
            self.inputs = qp.StlsInput(coupling, degeneracy, theory)
        else:
            self.throwInputFormatError("StlsInput:")
        self.isDefault = True
        
    # Set initial guess for chemical potential
    def setChemicalPotentialGuess(self, num1, num2):
        goodInput = ( isinstance(num1, (int, float)) and \
                      isinstance(num2, (int, float)) )
        if (goodInput) :
            self.inputs.setChemicalPotentialGuess(num1, num2)
            self.isDefault = False
        else:
            self.throwInputFormatError("setChemicalPotentialGuess")

    # Set miminum error for convergence
    def setErrMin(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setErrMin(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setErrMin")

    # Set iet mapping
    def setIETMapping(self, name):
        if (isinstance(name, str)):
            self.inputs.setIETMapping(name)
            self.isDefault = False
        else:
            self.throwInputFormatError("setIETMapping")
            
    # Set mixing parameter
    def setMixingParameter(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setMixingParameter(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setMixingParameter")

    # Set maximum number of iterations
    def setNIter(self, num):
        if (isinstance(num, int)):
            self.inputs.setNIter(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setNIter")
            
    # Set miminum error for convergence
    def setNMatsubara(self, num):
        if (isinstance(num, int)):
            self.inputs.setNMatsubara(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setNMatsubara")

    # Set output frequency
    def setOutIter(self, num):
        if (isinstance(num, int)):
            self.inputs.setOutIter(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setOutIter")

    # Set name for restart file
    def setRestartFileName(self, name):
        if (isinstance(name, str)):
            self.inputs.setRestartFileName(name)
            self.isDefault = False
        else:
            self.throwInputFormatError("setRestartFileName")

    # Set Wave-vector grid resolution
    def setWaveVectorGridRes(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setWaveVectorGridRes(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setWaveVectorGridRes")

    # Set Wave-vector grid cutoff
    def setWaveVectorGridCutoff(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setWaveVectorGridCutoff(num)
            self.isDefault = False
        else:
            self.throwInputFormatError("setWaveVectorGridCutoff")
                    
    # Print content of input structure
    def print(self):
        self.inputs.print()


class QstlsInput(StlsInput) :

    # Constructor
    def __init__(self, coupling, degeneracy, theory):
        goodInput = ( isinstance(coupling, (int, float)) and \
                      isinstance(degeneracy, (int, float)) and \
                      isinstance(theory, str) )
        if (goodInput):
            self.inputs = qp.StlsInput(coupling, degeneracy, theory)
            self.qinputs = qp.QstlsInput()
        else:
            self.throwInputFormatError("QstlsInput:")
        self.isDefault = True
    
    # Set name for file storing the fixed component of the auxiliary density
    # response
    def setFixedAdrFileName(self, name):
        if (isinstance(name, str)):
            self.qinputs.setFixedFileName(name)
            self.isDefault = False
        else:
            self.throwInputFormatError("setFixedAdrFileName")

    # Print content of input structure
    def print(self):
        self.inputs.print()
        self.qinputs.print()
