#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#----------------------------------------------------------------------------
#         PROJECT S-TIMATOR
#
# Unit test for S-timator Parser class
# Antonio Ferreira July 2006
#----------------------------------------------------------------------------

"""Unit test for modelparser.py"""

import sys
import os.path
import math

#append parent directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
import modelparser
import unittest

valid_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]

@ 3.4 pi = 2*pi

genomesize = 50 #should be enough
generations = 400

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt

"""

invalidIdDef_modelText ="""
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square
pipipi = pois  #this is an error'

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]

@ 3.4 pi = 2*pi

genomesize = 50 #should be enough
generations = 400

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt

"""

invalidIdFind_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1 + kkk ]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]

@ 3.4 pi = 2*pi

genomesize = 50 #should be enough
generations = 400

timecourse my file.txt  # this is a timecourse filename
timecourse anotherfile.txt
#timecourse stillanotherfile.txt

"""

invalidIdRate_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]

"""

invalidSyntaxComment_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG
But this is an invalid line

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2
pi   = 3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]
"""

invalidSyntaxDef_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2
pi   = 2*3.1416)
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]
"""

invalidSyntaxRate_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2)

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2
pi   = 2*3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]
"""

repeatedDeclaration_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2
reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = 2*3.1416
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]
"""

overflow_modelText = """
#This is an example of a valid model:

variables: SDLTSH TSH2 MG

Glx1 : TSH2  + MG -> SDLTSH, rate = Vmax1*TSH2*MG / ((KmMG+MG2)*(KmTSH2+TSH2))

reaction Glx2 : SDLTSH ->  , rate = Vmax2*SDLTSH / (Km2 + SDLTSH) #reaction 2

pi   = pi*1e100**10000
pi2  = 2*pi
pipi = pi**2  #this is pi square

find Vmax1 in [1e-9, 1e-3]
find   KmMG  in [1e-5, 1]
find KmTSH2 in [1e-5, pi/pi]

find Km2   in [1e-5, 1]
find Vmax2 in [1e-9, 1e-3]
"""

class LegalInputTests(unittest.TestCase):                          
    def setUp(self):
        self.parser = modelparser.StimatorParser()
        self.textlines = valid_modelText.split("\n")

    def testLegalModel(self):                          
        """parse() should give no error with legal model input"""
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, None)
        self.assertEqual(self.parser.variables, ['SDLTSH', 'TSH2', 'MG'])
        self.assertEqual(len(self.parser.constants), 3)
        self.assertEqual(self.parser.constants['pi'], 3.1416)
        self.assertEqual(self.parser.constants['pi2'], 2*3.1416)
        self.assertEqual(self.parser.constants['pipi'], 3.1416*3.1416)
        self.assertEqual(len(self.parser.parameters), 5)
        self.assertEqual(self.parser.parameters[0], ('Vmax1', 1e-9, 1e-3))
        self.assertEqual(self.parser.parameters[1], ('KmMG', 1e-5, 1.0))
        self.assertEqual(self.parser.parameters[2], ('KmTSH2', 1e-5, 1.0))
        self.assertEqual(self.parser.parameters[3], ('Km2', 1e-5, 1.0))
        self.assertEqual(self.parser.parameters[4], ('Vmax2', 1e-9, 1e-3))
        self.assertEqual(self.parser.timecourses, ['my file.txt', 'anotherfile.txt'])
        self.assertEqual(len(self.parser.rates), 2)
        self.assertEqual(self.parser.rates[0]['name'], 'Glx1')
        self.assertEqual(self.parser.rates[0]['irreversible'], True)
        self.assertEqual(self.parser.rates[0]['reagents'], [('TSH2', 1.0), ('MG', 1.0)])
        self.assertEqual(self.parser.rates[0]['products'], [('SDLTSH', 1.0)])
        self.assertEqual(self.parser.rates[0]['rate'], 'Vmax1*TSH2*MG / ((KmMG+MG)*(KmTSH2+TSH2))')
        self.assertEqual(self.parser.rates[1]['name'], 'Glx2')
        self.assertEqual(self.parser.rates[1]['irreversible'], True)
        self.assertEqual(self.parser.rates[1]['reagents'], [('SDLTSH', 1.0)])
        self.assertEqual(self.parser.rates[1]['products'], [])
        self.assertEqual(self.parser.rates[1]['rate'], 'Vmax2*SDLTSH / (Km2 + SDLTSH)')
        self.assertEqual(len(self.parser.atdefs), 1)
        self.assertEqual(self.parser.atdefs[0], (3.4, 'pi', 2*3.1416))
        self.assertEqual(len(self.parser.stoichmatrixrows), 3)
        self.assertEqual(len(self.parser.stoichmatrixrows[0]), 2)
        self.assertEqual(len(self.parser.stoichmatrixrows[1]), 1)
        self.assertEqual(len(self.parser.stoichmatrixrows[2]), 1)
        self.assertEqual(self.parser.stoichmatrixrows[0]['Glx1'], 1.0)
        self.assertEqual(self.parser.stoichmatrixrows[0]['Glx2'], -1.0)
        self.assertEqual(self.parser.stoichmatrixrows[1]['Glx1'], -1.0)
        self.assertEqual(self.parser.stoichmatrixrows[2]['Glx1'], -1.0)
        #modelparser.printParserResults(self.parser) #should this be called in an unit test?

class IllegalInputTests(unittest.TestCase):                          
    def setUp(self):
        self.parser = modelparser.StimatorParser()

    def testUnknownIdinDef(self):                          
        """parse() should give an error with unknown identifier in a constant definition"""
        self.textlines = invalidIdDef_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, "NameError : name 'pois' is not defined")
        self.assertEqual(self.parser.errorline, 12)
        self.assertEqual(self.parser.errorstart, 9)
        self.assertEqual(self.parser.errorend, 13)
        
    def testUnknownIdinFind(self):                          
        """parse() should give an error with unknown identifier in a parameter declaration"""
        self.textlines = invalidIdFind_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, "NameError : name 'kkk' is not defined")

    def testUnknownIdinRate(self):                          
        """parse() should give an error with unknown identifier in a rate definition"""
        self.textlines = invalidIdRate_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, "NameError : name 'MG2' is not defined")

    def testInvalidSyntaxComment(self):                          
        """parse() should give an error with a comment with missing #"""
        self.textlines = invalidSyntaxComment_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, "Invalid syntax")

    def testInvalidSyntaxDefExpression(self):                          
        """parse() should give an error with a constant definition without closing )"""
        self.textlines = invalidSyntaxDef_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assert_(self.parser.error.startswith('SyntaxError'))

    def testInvalidSyntaxRateExpression(self):                          
        """parse() should give an error with a rate definition without closing )"""
        self.textlines = invalidSyntaxRate_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assert_(self.parser.error.startswith('SyntaxError'))

    def testOverflow(self):                          
        """parse() should catch overflow errors"""
        self.textlines = overflow_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assert_(self.parser.error.startswith('OverflowError'))

    def testRepeatedDeclaration(self):                          
        """parse() should return an error if a declaration is repeated"""
        self.textlines = repeatedDeclaration_modelText.split("\n")
        self.parser.parse(self.textlines)
        self.assertEqual(self.parser.error, 'Repeated declaration')
        
if __name__ == "__main__":
    unittest.main()   

