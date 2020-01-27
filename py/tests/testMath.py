# -*- coding: utf-8 -*-
# This is the test of all C++ lib/high-precision/MathFunctions.hpp exported to python via py/high-precision/_math.cpp
# (C) 2015 Anton Gladky <gladk@debian.org>
# (C) 2019 Janek Kozicki

import unittest, math, sys
import yade
from yade import math as mth

import testMathHelper as mpmath
from   testMathHelper import mpc

class SimpleTests(unittest.TestCase):
	def setUp(self):
		self.printedAlready=set()
		self.extraName=""
		# the non-boost MPFR is not used in yade, that was just for more broader testing in minieigen-real, it used self.extraName="nb"
		#if(str("${LIBTOTEST}")[-2:] == "nb"):   extraName="nb" # non-boost MPFR
		self.nonBoostMPFR=False
		#if(str("${LIBTOTEST}")[-4:] == "BBFL"): extraName="_b" # boost cpp_bin_float
		if(yade.config.highPrecisionName.startswith("PrecisionBoost")):
			self.extraName="_b"
		self.digs0=yade.config.highPrecisionDecimalPlaces
		mpmath.mp.dps=self.digs0+1
		# tolerance = 1.2×10⁻ᵈ⁺¹, where ᵈ==self.digs0
		# so basically we store one more decimal digit, and expect one less decimal digit. That amounts to ignoring one (actually two) least significant digits.
		self.tolerance=(mpmath.mpf(10)**(-self.digs0+1))*mpmath.mpf("1.2")
		self.bits=mpmath.ceil(mpmath.mpf(self.digs0)/(mpmath.log(2)/mpmath.log(10)))+1
		# mpmath has 5 more internal bits
		self.expectedEpsilon=(2**5)*mpmath.eps()
		if(self.digs0 == 6): # float case
			self.bits=24
			self.expectedEpsilon=1.1920928955078125e-07
		if(self.digs0 == 15): # double case
			self.bits=53
			self.expectedEpsilon=2.220446049250313e-16
		if(self.digs0 == 18): # long double case
			self.bits=64
			self.expectedEpsilon=mpmath.mpf('1.084202172485504433993e-19')
		if(self.digs0 == 33): # float128 case
			self.bits=113
			self.expectedEpsilon=mpmath.mpf('1.925929944387235853055977942584926994e-34')
		if(yade.config.highPrecisionMpmath):
			self.maxval=(mpmath.mpf(1)-self.expectedEpsilon)*mpmath.power(2,mth.max_exp2)
		else:
			import sys
			self.maxval=sys.float_info.max

		# If failures appear and function is not broken then increase tolerance a little.
		self.defaultTolerances={
			#  function decimal places : tolerance factor. Each "10" corresponds to single wrong decimal place.
			#
			#                 float   double    long double float128        MPFR_100        MPFR_150     cpp_bin_float_100  cpp_bin_float_150
			#
			   "acos"      : {"6":100 , "15":100 , "18":100  , "33":1000   , "100":1000  , "150" :1000  , "100_b" :1000    , "150_b" :1000   }
			 , "atanh"     : {"6":100 , "15":100 , "18":100  , "33":1000   , "100":1000  , "150" :1000  , "100_b" :1000    , "150_b" :1000   }
			 , "acosh"     : {"6":100 , "15":100 , "18":100  , "33":1000   , "100":1000  , "150" :1000  , "100_b" :1000    , "150_b" :1000   }
			 , "asin"      : {"6":100 , "15":100 , "18":100  , "33":1000   , "100":1000  , "150" :1000  , "100_b" :1000    , "150_b" :1000   }

			# maybe the error lies in  mpmath, b ecause everything is compared with it.
			 , "sin"       : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "cos"       : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "tan"       : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }

			 , "csin"      : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "ccos"      : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "ctan"      : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }

			 , "cexp"      : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }
			 , "exp"       : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }
			 , "exp2"      : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }
			 , "expm1"     : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }
			 , "cosh"      : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }
			 , "sinh"      : {"6":10  , "15":10  , "18":10   , "33":10     , "100":10    , "150" :10    , "100_b" :100     , "150_b" :100    }

			 , "ccosh"     : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "csinh"     : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }
			 , "ctanh"     : {"6":100 , "15":100 , "18":20000, "33":4000   , "100":80000 , "150" :80000 , "100_b" :800000  , "150_b" :800000 }

			 , "clog"      : {"6":1000, "15":1000, "18":100  , "33":100    , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }
			 , "log"       : {"6":1000, "15":1000, "18":100  , "33":100    , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }
			 , "log10"     : {"6":1000, "15":1000, "18":100  , "33":100    , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }
			 , "log1p"     : {"6":1000, "15":1000, "18":100  , "33":100    , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }
			 , "log2"      : {"6":1000, "15":1000, "18":100  , "33":100    , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }

			 , "pow"       : {"6":5   , "15":10  , "18":50   , "33":50     , "100":50    , "150" :50    , "100_b" :50      , "150_b" :50     }
			 , "sqrt"      : {"6":5   , "15":10  , "18":50   , "33":50     , "100":50    , "150" :50    , "100_b" :50      , "150_b" :50     }

			 , "lgamma"    : {"6":100 , "15":100 , "18":1000 , "33":10000  , "100":100000, "150" :100000, "100_b" :1000000 , "150_b" :1000000}
			 , "tgamma"    : {"6":100 , "15":100 , "18":1000 , "33":10000  , "100":100000, "150" :100000, "100_b" :1000000 , "150_b" :1000000}
			 , "erfc"      : {"6":100 , "15":100 , "18":2000 , "33":20000  , "100":200000, "150" :200000, "100_b" :4000000 , "150_b" :8000000}
			 , "erf"       : {"6":100 , "15":100 , "18":5    , "33":5      , "100":100   , "150" :100   , "100_b" :100     , "150_b" :100    }

			 , "modf"      : {"6":10  , "15":100 , "18":5000 , "33":300000 , "100":10000 , "150" :100000, "100_b" :10000   , "150_b" :10000  }
			 , "fmod"      : {"6":10  , "15":100 , "18":5000 , "33":10000  , "100":10000 , "150" :100000, "100_b" :10000   , "150_b" :10000  }
			 , "remainder" : {"6":100 , "15":5000, "18":5000 , "33":10000  , "100":10000 , "150" :100000, "100_b" :10000   , "150_b" :10000  }
			 , "remquo"    : {"6":100 , "15":5000, "18":5000 , "33":10000  , "100":10000 , "150" :100000, "100_b" :10000   , "150_b" :10000  }
			 , "fma"       : {"6":10  , "15":100 , "18":10   , "33":10     , "100":100   , "150" :100   , "100_b" :100     , "150_b" :1000   }
			 }

	def getDefaultTolerance(self,name):
		dictForThisFunc = self.defaultTolerances[name]
		key = str(self.digs0)+self.extraName
		if(key in dictForThisFunc):
			return dictForThisFunc[key]*self.tolerance
		self.assertGreaterEqual(self.digs0 , 39) ## lower than 39 digits are all hardware precision: 6, 15, 18, 33 digits.
		low = dictForThisFunc["100"+self.extraName]
		high= dictForThisFunc["150"+self.extraName]
		import numpy
		return numpy.interp(self.digs0,[100,150],[low,high])

	def printOnce(self,functionName,a):
		if(functionName and (functionName not in self.printedAlready) and (not mpmath.isnan(abs(a)))):
			self.printedAlready.add(functionName)
			print(functionName.ljust(15)+" : "+a.__repr__())

	def checkRelativeError(self,a,b,tol=None,functionName=None,isComplex=False):
		if(abs(b) <= self.maxval and abs(b) >= mth.smallest_positive()):
			#print("a= ",a," b= ",b," smallest=",mth.smallest_positive(), " maxval=",self.maxval)
			self.printOnce(functionName,a)
			if(mpmath.isnan(a)):
				if(functionName != "lgamma"): # lgamma triggers this warning too often.
					print("\033[93m Warning: \033[0m got NaN, cannot verify if: ",a," == " ,b, " that was for function: \033[93m ",functionName, " \033[0m")
			else:
				if(tol != None):
					if isComplex:
						self.assertLessEqual(abs( (mpmath.mpc(a)-mpmath.mpc(b))/mpmath.mpc(b) ),tol)
					else:
						self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/mpmath.mpf(b) ),tol)
				else:
					if(functionName in self.defaultTolerances):
						defaultToleranceForThisFunction = self.getDefaultTolerance(functionName)
						#print(defaultToleranceForThisFunction," ---- ",functionName)
						if isComplex:
							self.assertLessEqual(abs( (mpmath.mpc(a)-mpmath.mpc(b))/mpmath.mpc(b) ),defaultToleranceForThisFunction)
						else:
							self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/mpmath.mpf(b) ),defaultToleranceForThisFunction)
					else:
						if isComplex:
							self.assertLessEqual(abs( (mpmath.mpc(a)-mpmath.mpc(b))/mpmath.mpc(b) ),self.tolerance)
						else:
							self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/mpmath.mpf(b) ),self.tolerance)
		else:
			print("Skipping ",functionName," check, the builtin number: ", a, " cannot have value outside of its possible repesentation: " , b, ", because it has only ",self.digs0," digits.")
	
	def checkRelativeComplexError(self,a,b,tol=None,functionName=None):
		self.printOnce(functionName,a)
		self.checkRelativeError(abs(a),abs(b),tol,functionName,True)

	def oneArgMathCheck(self,r):
		self.checkRelativeError(mth.sin(r),mpmath.sin(r),functionName="sin")
		self.checkRelativeError(mth.sinh(r),mpmath.sinh(r),functionName="sinh")
		self.checkRelativeError(mth.cos(r),mpmath.cos(r),functionName="cos")
		self.checkRelativeError(mth.cosh(r),mpmath.cosh(r),functionName="cosh")
		self.checkRelativeError(mth.tan(r),mpmath.tan(r),functionName="tan")
		self.checkRelativeError(mth.tanh(r),mpmath.tanh(r),functionName="tanh")
		# check math functions, but ensure that input arguments produce real (not complex) results
		self.checkRelativeError(mth.abs(r),abs(r),functionName="abs")
		self.checkRelativeError(mth.acos(r%1),mpmath.acos(r%1),functionName="acos")
		self.checkRelativeError(mth.acosh(abs(r)+1),mpmath.acosh(abs(r)+1),functionName="acosh")
		self.checkRelativeError(mth.asin(r%1),mpmath.asin(r%1),functionName="asin")
		self.checkRelativeError(mth.asinh(r),mpmath.asinh(r),functionName="asinh")
		self.checkRelativeError(mth.atan(r),mpmath.atan(r),functionName="atan")
		self.checkRelativeError(mth.atanh(r%1),mpmath.atanh(r%1),functionName="atanh")
		self.checkRelativeError(mth.cbrt(abs(r)),mpmath.cbrt(abs(r)),functionName="cbrt")
		self.assertEqual(mth.ceil(r),mpmath.ceil(r))
		self.checkRelativeError(mth.erf(r),mpmath.erf(r),functionName="erf")
		self.checkRelativeError(mth.erfc(r),mpmath.erfc(r),functionName="erfc")
		self.checkRelativeError(mth.exp(r),mpmath.exp(r),functionName="exp")
		self.checkRelativeError(mth.sqrt(abs(r)),mpmath.sqrt(abs(r)),functionName="sqrt")
		self.checkRelativeError(mth.exp2(r),mpmath.power(2,r),functionName="exp2")
		self.checkRelativeError(mth.expm1(r),mpmath.expm1(r),functionName="expm1")
		self.assertEqual(mth.floor(r),mpmath.floor(r))
		#print(mth.ilogb(r).__repr__()) # ilogb is not present in mpmath
		self.checkRelativeError(mth.lgamma(r),mpmath.log(abs(mpmath.gamma(r))),functionName="lgamma")
		self.checkRelativeError(mth.log(abs(r)+self.tolerance),mpmath.log(abs(r)+self.tolerance),functionName="log")
		self.checkRelativeError(mth.log10(abs(r)+self.tolerance),mpmath.log10(abs(r)+self.tolerance),functionName="log10")
		self.checkRelativeError(mth.log1p(abs(r)+self.tolerance),mpmath.log(1+abs(r)+self.tolerance),functionName="log1p")
		self.checkRelativeError(mth.log2(abs(r)+self.tolerance),mpmath.log(abs(r)+self.tolerance)/mpmath.log(2),functionName="log2")
		#print(mth.logb(r).__repr__()) # logb is not present in mpmath
		self.assertEqual(mth.rint(r),round(r))
		self.assertTrue((mth.round(r)==round(r)) or (r%1==0.5)) # ignore rounding 0.5 up or down.
		self.checkRelativeError(mth.tgamma(r),mpmath.gamma(r),functionName="tgamma")
		self.assertEqual(mth.trunc(abs(r)),int(abs(r)))

		self.checkRelativeError(mth.fabs(r),abs(r),functionName="fabs")

		pair = mth.frexp(abs(r))
		self.checkRelativeError(abs(r),pair[0]*mpmath.power(2,pair[1]),functionName="frexp")

		pair = mth.modf(abs(r))
		self.checkRelativeError(pair[0],(abs(r))%1,functionName="modf")
		self.assertEqual(pair[1],int(abs(r)))

		#self.assertEqual(mth.frexp(abs(r)),mth.frexp_c_test(abs(r)))
		#self.assertEqual(mth.modf(abs(r)),mth.modf_c_test(abs(r)))

		if(r==0):
			self.assertEqual(mth.sgn(r),0)
			self.assertEqual(mth.sign(r),0)
		if(r>0):
			self.assertEqual(mth.sgn(r),1)
			self.assertEqual(mth.sign(r),1)
		if(r<0):
			self.assertEqual(mth.sgn(r),-1)
			self.assertEqual(mth.sign(r),-1)

		self.checkCgalNumTraits(r)

	def checkCgalNumTraits(self,r):
		if(mth.testCgalNumTraits==False):
			print("Skipping test of CgalNumTraits")
			return
		self.assertEqual(mth.CGAL_Is_valid(r),True)
		if(r != 0):
			self.checkRelativeError(mth.CGAL_Square(r),mpmath.power(r,2),functionName="pow")
			self.checkRelativeError(mth.CGAL_Sqrt(abs(r)),mpmath.sqrt(abs(r)),functionName="sqrt")
			for kk in range(5):
				k=kk+1
				self.checkRelativeError(mth.CGAL_Kth_root(k,abs(r)),mpmath.power(abs(r),1/mpmath.mpf(k)),functionName="pow")
			# CGAL uses double for intervals
			interval = mth.CGAL_To_interval(r)
			self.checkRelativeError(r,interval[0],1e-14)
			self.checkRelativeError(r,interval[1],1e-14)
		self.assertEqual(mth.CGAL_Is_finite(r),True)
		if(r==0): self.assertEqual(mth.CGAL_Sgn(r),0)
		if(r> 0): self.assertEqual(mth.CGAL_Sgn(r),1)
		if(r< 0): self.assertEqual(mth.CGAL_Sgn(r),-1)
		self.assertEqual(mth.CGAL_Sgn(0),0)
		self.assertEqual(mth.CGAL_Sgn(2.5),1)
		self.assertEqual(mth.CGAL_Sgn(-2.3),-1)

	def testInfinityNaN(self):
		if(mth.hasInfinityNan==False):
			print("Skipping inf,nan regular test\n")
			print("\033[91m *** Warning: usually YADE needs Inf and NaN for most of the calculations. *** \033[0m")
			return
		self.assertEqual(mth.isinf   (mpmath.mpf( 1   )),False)
		self.assertEqual(mth.isinf   (mpmath.mpf('nan')),False)
		self.assertEqual(mth.isinf   (mpmath.mpf('inf')),True )
		self.assertEqual(mth.isnan   (mpmath.mpf( 1   )),False)
		self.assertEqual(mth.isnan   (mpmath.mpf('nan')),True )
		self.assertEqual(mth.isnan   (mpmath.mpf('inf')),False)
		self.assertEqual(mth.isfinite(mpmath.mpf( 1   )),True )
		self.assertEqual(mth.isfinite(mpmath.mpf('nan')),False)
		self.assertEqual(mth.isfinite(mpmath.mpf('inf')),False)

	def testCgalNumTraits(self):
		if(mth.testCgalNumTraits==False):
			print("Skipping test of CgalNumTraits")
			return
		self.checkCgalNumTraits(0)
		self.checkCgalNumTraits(0.5)
		self.checkCgalNumTraits(-1.5)
		self.checkCgalNumTraits(55.5)
		self.assertEqual(mth.CGAL_Is_valid (mpmath.mpf( 1   )),True )
		self.assertEqual(mth.CGAL_Is_valid (mpmath.mpf('nan')),False)
		self.assertEqual(mth.CGAL_Is_valid (mpmath.mpf('inf')),True )
		self.assertEqual(mth.CGAL_Is_finite(mpmath.mpf( 1   )),True )
		self.assertEqual(mth.CGAL_Is_finite(mpmath.mpf('nan')),False)
		self.assertEqual(mth.CGAL_Is_finite(mpmath.mpf('inf')),False)
		self.assertEqual(mth.CGAL_simpleTest(),mpmath.mpf("3.0"))

	def twoArgMathCheck(self,r1,r2):
		self.checkRelativeComplexError(mth.sin (mpmath.mpc(r1,r2)),mpmath.sin (mpmath.mpc(r1,r2)),functionName="csin")
		self.checkRelativeComplexError(mth.sinh(mpmath.mpc(r1,r2)),mpmath.sinh(mpmath.mpc(r1,r2)),functionName="csinh")
		self.checkRelativeComplexError(mth.cos (mpmath.mpc(r1,r2)),mpmath.cos (mpmath.mpc(r1,r2)),functionName="ccos")
		self.checkRelativeComplexError(mth.cosh(mpmath.mpc(r1,r2)),mpmath.cosh(mpmath.mpc(r1,r2)),functionName="ccosh")
		self.checkRelativeComplexError(mth.tan (mpmath.mpc(r1,r2)),mpmath.tan (mpmath.mpc(r1,r2)),functionName="ctan")
		self.checkRelativeComplexError(mth.tanh(mpmath.mpc(r1,r2)),mpmath.tanh(mpmath.mpc(r1,r2)),functionName="ctanh")

		self.checkRelativeComplexError(mth.exp(mpmath.mpc(r1,r2)),mpmath.exp(mpmath.mpc(r1,r2)),functionName="cexp")
		self.checkRelativeComplexError(mth.log(mpmath.mpc(r1,r2)),mpmath.log(mpmath.mpc(r1,r2)),functionName="clog")

		self.checkRelativeComplexError(mth.abs  (mpmath.mpc(r1,r2)),abs(mpmath.mpc(r1,r2)),functionName="cabs")
		self.checkRelativeComplexError(mth.conj (mpmath.mpc(r1,r2)),mpmath.conj(mpmath.mpc(r1,r2)),functionName="cconj")
		self.checkRelativeComplexError(mth.real (mpmath.mpc(r1,r2)),r1,functionName="creal")
		self.checkRelativeComplexError(mth.imag (mpmath.mpc(r1,r2)),r2,functionName="cimag")

		self.checkRelativeError(mth.atan2(r1,r2),mpmath.atan2(r1,r2),functionName="atan2")
		self.checkRelativeError(mth.fmod(abs(r1),abs(r2)),mpmath.fmod(abs(r1),abs(r2)),functionName="fmod")
		self.checkRelativeError(mth.hypot(r1,r2),mpmath.hypot(r1,r2),functionName="hypot")
		self.checkRelativeError(mth.max(r1,r2),max(r1,r2),functionName="max")
		self.checkRelativeError(mth.min(r1,r2),min(r1,r2),functionName="min")
		self.checkRelativeError(mth.pow(abs(r1),r2),mpmath.power(abs(r1),r2),functionName="pow")
		self.checkRelativeError(mth.remainder(abs(r1),abs(r2)),abs(r1)-round(abs(r1)/abs(r2))*abs(r2),functionName="remainder")
		pair = mth.remquo(abs(r1),abs(r2))
		self.checkRelativeError(pair[0],abs(r1)-round(abs(r1)/abs(r2))*abs(r2),functionName="remquo")
		self.assertEqual(pair[1]%8, round(abs(r1/r2))%8)

		self.checkRelativeError(mth.ldexp(r1,int(r2)),mpmath.mpf(r1)*mpmath.power(2,int(r2)),functionName="ldexp")

	def threeArgMathCheck(self,r1,r2,r3):
		self.checkRelativeError(mth.fma(r1,r2,r3),(mpmath.mpf(r1)*r2)+r3,functionName="fma")

	def testMathFunctions(self):
		self.assertEqual(mth.defprec , self.bits )
		zz=mpmath.acos(0)
		#print(zz.__repr__())
		#print("zz:",hex(id(zz)))
		#print("mpmath:",hex(id(mpmath)))
		a=mth.Var()
		a.val=zz
		self.assertEqual(mpmath.mp.dps , self.digs0+1 )
		#print("---- a.val=",a.val.__repr__())
		#print("---- zz   =",zz   .__repr__())
		#print("---- DPS  =",mpmath.mp.dps)
		#print("---- abs  =",abs(mpmath.mpf(a.val-zz)))
		#print("---- 10** =",self.tolerance)
		self.checkRelativeError(a.val,zz)
		self.assertEqual(mth.IsInteger, 0 )
		self.assertEqual(mth.IsSigned, 1 )
		self.assertEqual(mth.IsComplex,  0)
		if(self.bits >= 64):
			self.assertEqual(mth.RequireInitialization, 1 )
		else:
			self.assertEqual(mth.RequireInitialization, 0 )
		self.assertGreaterEqual(mth.ReadCost, 1)
		self.assertGreaterEqual(mth.AddCost, 1)
		self.assertGreaterEqual(mth.MulCost, 1)
		self.checkRelativeError(mth.highest(),self.maxval,2.1)
		self.checkRelativeError(-mth.lowest(),self.maxval,2.1)
		self.checkRelativeError(mth.Pi(),mpmath.pi)
		self.checkRelativeError(mth.Euler(),mpmath.euler)
		self.checkRelativeError(mth.Log2(),mpmath.log(2))
		self.checkRelativeError(mth.Catalan(),mpmath.catalan)
		#print("mth.epsilon() ",mth.epsilon(),"  self.expectedEpsilon = ",self.expectedEpsilon)
		self.checkRelativeError(mth.epsilon(),self.expectedEpsilon,10)
		if(self.digs0 == 6): # exception for float
			self.assertLessEqual(mth.dummy_precision(),10e-6)
		else:
			self.checkRelativeError(mpmath.log(mth.dummy_precision()/mth.epsilon())/mpmath.log(10) , mpmath.mpf(self.digs0)/10 , 1.5 )
		for x in range(50):
			if(self.nonBoostMPFR): # this looks like a bug in /usr/include/eigen3/unsupported/Eigen/MPRealSupport !
				self.assertLessEqual(abs(mth.random()-0.5),0.5)
			else:
				self.assertLessEqual(abs(mth.random()    ),1.0)
		for aa in range(4):
			for bb in range(4):
				a = (aa-3)*5
				b = bb*10
				r = mth.random(a,a+b+1)
				r2= mth.random(a,a+b+1)
				r3= mth.random(a,a+b+5)
				#print("random=",r)
				self.assertLessEqual(r,a+b+1)
				self.assertGreaterEqual(r,a)
				self.assertFalse(mth.isMuchSmallerThan(r,1,mth.epsilon()))
				self.assertTrue(mth.isMuchSmallerThan(self.expectedEpsilon,1+abs(r),mth.epsilon()))
				self.assertTrue(mth.isEqualFuzzy(r+self.expectedEpsilon*0.01,r,mth.epsilon()))
				self.checkRelativeError(mth.toLongDouble(r),float(r), 1e-14) # FIXME - should be 1e-17, but python does not support that
				self.checkRelativeError(mth.toDouble(r),float(r), 1e-14)
				self.checkRelativeError(mth.toDouble(r),float(r), 1e-14)
				self.assertEqual(mth.toLong(r),int(r))
				self.assertEqual(mth.toInt(r),int(r))
				#
				#print(r.__repr__(),r2.__repr__(),r3.__repr__())
				self.oneArgMathCheck(r)
				self.oneArgMathCheck(r2)
				self.oneArgMathCheck(r3)
				self.twoArgMathCheck(r,r2)
				self.twoArgMathCheck(r,r3)
				self.twoArgMathCheck(r2,r3)
				self.threeArgMathCheck(r,r2,r3)

	def testArray(self):
		mth.testArray()

	def testBasicVariable(self):
		a=mth.Var()
		self.checkRelativeError(a.val,-71.23,0.01)
		a.val=10
		self.checkRelativeError(a.val,10)
		self.checkRelativeComplexError(a.cpl,-71.23+33.23j,0.01)
		a.cpl=mpmath.mpc("1","-1")
		self.checkRelativeComplexError(a.cpl,1-1j,1e-15)
		self.checkRelativeComplexError(a.cpl,mpmath.mpc("1","-1"))

	def thisTestsExceptionReal(self):
		a=mth.Var()
		a.val="13123-123123*123"

	def thisTestsExceptionComplex(self):
		a=mth.Var()
		a.cpl="13123-123123*123-50j"

	def testWrongInput(self):
		if(self.nonBoostMPFR): # this looks like another bug in /usr/include/mpreal.h
			print("skipping this test for non-boost /usr/include/mpreal.h")
			return
		# depending on backed Real use it throws TypeError or RuntimeError
		self.assertRaises(Exception,self.thisTestsExceptionReal)
		self.assertRaises(Exception,self.thisTestsExceptionComplex)

