#!/usr/bin/env python

import sys
import qupled.qupled as qp

class Input:

    # Constructor
    def __init__(self):
        self.qpInput = qp.Input()

    # Print content of input structure
    def print(self):
        self.qpInput.print()

    # Set coupling parameter
    def setCoupling(self, num):
        if (isinstance(num, (int, float))):
            self.qpInput.set("base.coupling", str(num))
        else:
            self.__throwInputFormatError("setCoupling")

    # Set degeneracy parameter
    def setDegeneracy(self, num):
        if (isinstance(num, (int, float))):
            self.qpInput.set("base.degeneracy", str(num))
        else:
            self.__throwInputFormatError("setDegeneracy")

    # Set minimum error for convergence
    def setConvergenceError(self, num):
        if (isinstance(num, (int, float))):
            self.qpInput.set("static.error", str(num))
        else:
            self.__throwInputFormatError("setConvergenceError")
            
    # Set mixing parameter
    def setMixing(self, num):
        if (isinstance(num, (int, float))):
            self.qpInput.set("static.mixing", str(num))
        else:
            self.__throwInputFormatError("setMixing")
        
    # Set number of OMP threads
    def setNumberOfThreads(self, num):
        if (isinstance(num, int)):
            self.qpInput.set("base.threads", str(num))
        else:
            self.__throwInputFormatError("setNumberOfThreads")
        
    # Set scheme for 2D integrals
    def setSchemeFor2DIntegrals(self, name):
        if (isinstance(num, string)):
            self.qpInput.set("base.int2DScheme", name)
        else:
            self.__throwInputFormatError("setSchemeFor2DIntegrals")
            
    # Set theory to solve
    def setTheory(self, name):
        if (isinstance(num, string)):
            self.qpInput.set("base.theory", name)
        else:
            self.__throwInputFormatError("setTheory")

    # Throw error messages
    def __throwInputFormatError(self, callerName):
        sys.exit(callerName + ": Wrong input format")


# static.waveVectorResolution = 0.1
# static.waveVectorCutoff = 10
# static.chemicalPotential = -10,10
# static.matsubara = 128
# static.iterations = 1000
# static.outputFrequency = 10
# stls.iet = standard
# stls.restart = 
# qstls.useStatic = 0
# qstls.fixed = 

