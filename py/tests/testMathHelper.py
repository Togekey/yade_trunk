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
	def __init__(self):
		pass
	class mpf:
		def __init__(self,realpart):
			self.r = float(realpart)
		def __float__(self):
			return float(self.r)
		def __pow__(self,p):
			return float(self.r**float(p))
		def __truediv__(self,p):
			return float(self.r/float(p))
		def __rtruediv__(self,p):
			return float(float(p)/self.r)
		def __rmul__(self,p):
			return float(self.r*float(p))
		def __sub__(self,b):
			return float(self.r - float(b))
	class mpc:
		def __init__(self,realpart, imagpart=None):
			if(imagpart.__class__ == str):
				self.i = float(imagpart)
			elif(imagpart == None):
				self.i = 0
			else:
				raise TypeError
			if(realpart.__class__ == complex or realpart.__class__ == SuperMinimalMath.mpc):
				self.r = complex(realpart).real
				self.i = complex(realpart).imag
			elif(realpart.__class__ == float or realpart.__class__ == int or realpart.__class__ == SuperMinimalMath.mpf or realpart.__class__ == str):
				self.r = float(realpart)
			else:
				raise TypeError
		def __complex__(self):
			return complex(self.r + self.i*1j)
		def __sub__(self,b):
			return complex(self) - complex(b)
		def __truediv__(self,p):
			return complex(self)/complex(p)
		def __rtruediv__(self,p):
			return complex(p)/complex(self)
		def __abs__(self):
			return abs(complex(self))
