# encoding: utf-8
##########################################################################
#  2019      Janek Kozicki                                               #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
"""
This module is dedicated to testing python ↔ C++ interface and exposing
to ``yade --test`` stuff which isn't exported to python in usual ways.
It covers:

#. The testing of python / C++ interface, e.g. virtual functions inheritance etc.
#. Handling of numerical precision and numerical conversion for python ↔ C++
#. The pure C++ functions, e.g. math or tensor functions not exposed to python
#. Allows to call internal C++ testing of functions not exposed to python

Example session:

.. ipython::

   In [1]: from yade.testCppPy import *

   In [1]: testPiPython()

   In [1]: testPiFloat()

   In [1]: testPiDouble()

   In [1]: testPiLongDouble()

   In [1]: res=[ testPiPython() , testPiFloat() , testPiDouble() , testPiLongDouble() ]

   In [1]: print("\\n".join(["%.16f " % r[1] for r in res]))

"""
# all C++ functions are accessible now:
from yade._testCppPy import *

def testPiPython():
	"""
	This function returns a string and python numerical representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

	:return: ``(string(pi),pi)`` from python math.pi
	"""
	import math
	print("python::testPiPython()                                           returns pi = %.16f " % math.pi)
	return ( "%.16f" % math.pi , math.pi )

