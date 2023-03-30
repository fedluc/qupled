#!/usr/bin/env python

import sys
import qupled.qupled as qp

class Base:

    # Constructor
    def __init__(self, coupling, degeneracy, theory):
        goodInput = ( isinstance(coupling, (int, float)) and \
                      isinstance(degeneracy, (int, float)) and \
                      isinstance(theory, str) )
        if (goodInput):
            self.inputs = qp.Input(coupling, degeneracy, theory)
        else:
            self.throwInputFormatError("Input:")

    # Set coupling parameter
    def setCoupling(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setCoupling(num)
        else:
            self.throwInputFormatError("setCoupling")

    # Set degeneracy parameter
    def setDegeneracy(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setDegeneracy(num)
        else:
            self.throwInputFormatError("setDegeneracy")

    # Set theory to solve
    def setTheory(self, name):
        if (isinstance(name, str)):
            self.inputs.setTheory(name)
        else:
            self.throwInputFormatError("setTheory")
            
    # Set number of OMP threads
    def setNThreads(self, num):
        if (isinstance(num, int)):
            self.inputs.setNThreads(num)
        else:
            self.throwInputFormatError("setNThreads")
    
    # Set scheme for 2D integrals
    def setInt2DScheme(self, name):
        if (isinstance(name, str)):
            self.inputs.setInt2DScheme(name)
        else:
            self.throwInputFormatError("setInt2DScheme")

    # Get coupling parameter
    def getCoupling(self):
        return self.inputs.getCoupling()

    # Get degeneracy parameter
    def getDegeneracy(self):
        return self.inputs.getDegeneracy()

    # Get theory to be solved
    def getTheory(self):
        return self.inputs.getTheory()
    
    # Get number of OMP threads
    def getNThreads(self):
        return self.inputs.getNThreads()
    
    # Get scheme for 2D integrals
    def getInt2DScheme(self):
        return self.inputs.getInt2DScheme()

    # Print content of input structure
    def print(self):
        self.inputs.print()

    # Compare to another Input object
    def isEqual(self, obj):
        return self.inputs.isEqual(obj.inputs)

    # Throw error messages
    def throwInputFormatError(self, callerName):
        sys.exit(callerName + ": Wrong input format")


class Stls(Base):

    # Constructor
    def __init__(self, coupling, degeneracy, theory):
        goodInput = ( isinstance(coupling, (int, float)) and \
                      isinstance(degeneracy, (int, float)) and \
                      isinstance(theory, str) )
        if (goodInput):
            self.inputs = qp.StlsInput(coupling, degeneracy, theory)
        else:
            self.throwInputFormatError("StlsInput:")
        
    # Set initial guess for chemical potential
    def setChemicalPotentialGuess(self, num1, num2):
        goodInput = ( isinstance(num1, (int, float)) and \
                      isinstance(num2, (int, float)) )
        if (goodInput) :
            self.inputs.setChemicalPotentialGuess(num1, num2)
        else:
            self.throwInputFormatError("setChemicalPotentialGuess")

    # Set miminum error for convergence
    def setErrMin(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setErrMin(num)
        else:
            self.throwInputFormatError("setErrMin")

    # Set iet mapping
    def setIETMapping(self, name):
        if (isinstance(name, str)):
            self.inputs.setIETMapping(name)
        else:
            self.throwInputFormatError("setIETMapping")
            
    # Set mixing parameter
    def setMixingParameter(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setMixingParameter(num)
        else:
            self.throwInputFormatError("setMixingParameter")

    # Set maximum number of iterations
    def setNIter(self, num):
        if (isinstance(num, int)):
            self.inputs.setNIter(num)
        else:
            self.throwInputFormatError("setNIter")
            
    # Set miminum error for convergence
    def setNMatsubara(self, num):
        if (isinstance(num, int)):
            self.inputs.setNMatsubara(num)
        else:
            self.throwInputFormatError("setNMatsubara")

    # Set output frequency
    def setOutIter(self, num):
        if (isinstance(num, int)):
            self.inputs.setOutIter(num)
        else:
            self.throwInputFormatError("setOutIter")

    # Set name for restart file
    def setRestartFileName(self, name):
        if (isinstance(name, str)):
            self.inputs.setRestartFileName(name)
        else:
            self.throwInputFormatError("setRestartFileName")

    # Set Wave-vector grid resolution
    def setWaveVectorGridRes(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setWaveVectorGridRes(num)
        else:
            self.throwInputFormatError("setWaveVectorGridRes")

    # Set Wave-vector grid cutoff
    def setWaveVectorGridCutoff(self, num):
        if (isinstance(num, (int, float))):
            self.inputs.setWaveVectorGridCutoff(num)
        else:
            self.throwInputFormatError("setWaveVectorGridCutoff")

    # Get initial guess for chemical potential
    def getChemicalPotentialGuess(self):
        return self.inputs.getChemicalPotentialGuess()

    # Get miminum error for convergence
    def getErrMin(self):
        return self.inputs.getErrMin()

    # Get iet mapping
    def getIETMapping(self):
        return self.inputs.getIETMapping()
            
    # Get mixing parameter
    def getMixingParameter(self):
        return self.inputs.getMixingParameter()

    # Get maximum number of iterations
    def getNIter(self):
        return self.inputs.getNIter()
            
    # Get miminum error for convergence
    def getNMatsubara(self):
        return self.inputs.getNMatsubara()

    # Get output frequency
    def getOutIter(self):
        return self.inputs.getOutIter()

    # Get name for restart file
    def getRestartFileName(self):
        return self.inputs.getRestartFileName()

    # Get Wave-vector grid resolution
    def getWaveVectorGridRes(self):
        return self.inputs.getWaveVectorGridRes()

    # Get Wave-vector grid cutoff
    def getWaveVectorGridCutoff(self):
        return self.inputs.getWaveVectorGridCutoff()

    # Compare to another StlsInput object
    def isEqual(self, obj):
        return self.inputs.isEqual(obj.inputs)
    
    # Print content of input structure
    def print(self):
        self.inputs.print()


class Qstls(Stls) :

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
    
    # Set name for file storing the fixed component of the auxiliary density
    # response
    def setFixedAdrFileName(self, name):
        if (isinstance(name, str)):
            self.qinputs.setFixedFileName(name)
        else:
            self.throwInputFormatError("setFixedAdrFileName")

    # Get name for file storing the fixed component of the auxiliary density
    # response
    def getFixedAdrFileName(self):
        return self.qinputs.setFixedFileName()

    # Compare to another QstlsInput object
    def isEqual(self, obj):
        return ( self.inputs.isEqual(obj.inputs) and
                 self.qinputs.isEqual(obj.qinputs) )

    # Print content of input structure
    def print(self):
        self.inputs.print()
        self.qinputs.print()
