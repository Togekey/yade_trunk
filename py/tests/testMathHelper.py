#!/usr/bin/python
# -*- coding: utf-8 -*-
##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

# When mpmath is not required implement a super-minimal version of mpmath, so that the tests will work and pretend that they use mpmath
# So this file is in fact to avoid python3-mpmath dependency when high-precision is not used.

class MP:
	dps=None
class SuperMinimalMath:
	mp=MP
	mpf=float
	# mpc needs to accept two string arguments, apart from that it's just a complex number
	class mpc(complex):
		def __new__(cls, *args, **kwargs):
			super_new = super(SuperMinimalMath.mpc, cls).__new__
			if super_new is object.__new__:
				return super_new(cls)
			if len(args)==2:
				return super_new(cls, float(args[0]), float(args[1]))
			return super_new(cls, *args, **kwargs)

