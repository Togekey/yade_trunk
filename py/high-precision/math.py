# encoding: utf-8
##########################################################################
#  2019      Janek Kozicki                                               #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
"""
This module exports all C++ math functions available in yade.
"""
# all C++ functions are accessible now:
from yade._math import *

def testPiPython():
	"""
	This function returns a string and python numerical representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

	:return: ``(string(pi),pi)`` from python math.pi
	"""
	import math
	return ( "%.16f" % math.pi , math.pi )

