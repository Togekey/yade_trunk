# encoding: utf-8
##########################################################################
#  2019      Janek Kozicki                                               #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
"""
This python module exposes all C++ math functions for Real and Complex type to python.
In fact it sort of duplicates ``import math``, ``import cmath`` or ``import mpmath``.
This module has four purposes:

1. to reliably test all C++ math functions of arbitrary precision Real and Complex types against mpmath.
2. to test Eigen NumTraits
3. to test CGAL NumTraits
4. To allow writing python math code in a way that mirrors C++ math code in yade. As a bonus it
   will be faster than mpmath because mpmath is a purely python library (which was one of the main
   difficulties when writing :ysrc:`lib/high-precision/ToFromPythonConverter.hpp`)

If another math function is needed it should be added to following files:

1. :ysrc:`lib/high-precision/MathFunctions.hpp`
2. :ysrc:`py/high-precision/_math.cpp`
3. :ysrc:`py/tests/testMath.py`

"""
# all C++ functions are accessible now:
from yade._math import *

