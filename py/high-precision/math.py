# encoding: utf-8
##########################################################################
#  2019      Janek Kozicki                                               #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
"""
This python module exposes all C++ math functions for Real and Complex type to python.
In fact it sort of duplicates 'import math', 'import cmath' or 'import mpmath'.
This module has three purposes:

1. to reliably test all C++ math functions of arbitrary precision Real and Complex types against mpmath.
2. to test Eigen NumTraits
3. To allow writing python math code in a way that mirrors C++ math code in yade. As a bonus it
   will be faster than mpmath because mpmath is a purely python library (which was one of the main
   difficulties when writing lib/high-precision/ToFromPythonConverter.hpp)

"""
# all C++ functions are accessible now:
from yade._math import *

