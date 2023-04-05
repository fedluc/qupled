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

    # Compare to another Stls object
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

    # Compare to another Qstls object
    def isEqual(self, obj):
        return ( self.inputs.isEqual(obj.inputs) and
                 self.qinputs.isEqual(obj.qinputs) )

    # Print content of input structure
    def print(self):
        self.inputs.print()
        self.qinputs.print()
