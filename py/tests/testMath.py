# -*- coding: utf-8 -*-
# This is the test of all C++ lib/high-precision/MathFunctions.hpp exported to python via py/high-precision/_math.cpp
# (C) 2015 Anton Gladky <gladk@debian.org>
# (C) 2019 Janek Kozicki

import unittest, math, sys
import yade
from yade import math as mne

import testMathHelper as mpmath
from   testMathHelper import mpc

class SimpleTests(unittest.TestCase):
	def setUp(self):
		self.printCount=0
		self.everyNth=100
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
			self.maxval=(mpmath.mpf(1)-self.expectedEpsilon)*mpmath.power(2,mne.max_exp2)
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

	def checkRelativeError(self,a,b,tol=None,functionName=None,isComplex=False):
		if(abs(b) <= self.maxval and abs(b) >= mne.smallest_positive()):
			#print("a= ",a," b= ",b," smallest=",mne.smallest_positive(), " maxval=",self.maxval)
			self.printCount+=1
			if(functionName and (self.printCount % self.everyNth == 0)):
				print(functionName.ljust(15)+" : "+a.__repr__())
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
		self.checkRelativeError(abs(a),abs(b),tol,functionName,True)

	def oneArgMathCheck(self,r):
		self.checkRelativeError(mne.sin(r),mpmath.sin(r),functionName="sin")
		self.checkRelativeError(mne.sinh(r),mpmath.sinh(r),functionName="sinh")
		self.checkRelativeError(mne.cos(r),mpmath.cos(r),functionName="cos")
		self.checkRelativeError(mne.cosh(r),mpmath.cosh(r),functionName="cosh")
		self.checkRelativeError(mne.tan(r),mpmath.tan(r),functionName="tan")
		self.checkRelativeError(mne.tanh(r),mpmath.tanh(r),functionName="tanh")
		# check math functions, but ensure that input arguments produce real (not complex) results
		self.checkRelativeError(mne.abs(r),abs(r),functionName="abs")
		self.checkRelativeError(mne.acos(r%1),mpmath.acos(r%1),functionName="acos")
		self.checkRelativeError(mne.acosh(abs(r)+1),mpmath.acosh(abs(r)+1),functionName="acosh")
		self.checkRelativeError(mne.asin(r%1),mpmath.asin(r%1),functionName="asin")
		self.checkRelativeError(mne.asinh(r),mpmath.asinh(r),functionName="asinh")
		self.checkRelativeError(mne.atan(r),mpmath.atan(r),functionName="atan")
		self.checkRelativeError(mne.atanh(r%1),mpmath.atanh(r%1),functionName="atanh")
		self.checkRelativeError(mne.cbrt(abs(r)),mpmath.cbrt(abs(r)),functionName="cbrt")
		self.assertEqual(mne.ceil(r),mpmath.ceil(r))
		self.checkRelativeError(mne.erf(r),mpmath.erf(r),functionName="erf")
		self.checkRelativeError(mne.erfc(r),mpmath.erfc(r),functionName="erfc")
		self.checkRelativeError(mne.exp(r),mpmath.exp(r),functionName="exp")
		self.checkRelativeError(mne.sqrt(abs(r)),mpmath.sqrt(abs(r)),functionName="sqrt")
		self.checkRelativeError(mne.exp2(r),mpmath.power(2,r),functionName="exp2")
		self.checkRelativeError(mne.expm1(r),mpmath.expm1(r),functionName="expm1")
		self.assertEqual(mne.floor(r),mpmath.floor(r))
		#print(mne.ilogb(r).__repr__()) # ilogb is not present in mpmath
		self.checkRelativeError(mne.lgamma(r),mpmath.log(abs(mpmath.gamma(r))),functionName="lgamma")
		self.checkRelativeError(mne.log(abs(r)+self.tolerance),mpmath.log(abs(r)+self.tolerance),functionName="log")
		self.checkRelativeError(mne.log10(abs(r)+self.tolerance),mpmath.log10(abs(r)+self.tolerance),functionName="log10")
		self.checkRelativeError(mne.log1p(abs(r)+self.tolerance),mpmath.log(1+abs(r)+self.tolerance),functionName="log1p")
		self.checkRelativeError(mne.log2(abs(r)+self.tolerance),mpmath.log(abs(r)+self.tolerance)/mpmath.log(2),functionName="log2")
		#print(mne.logb(r).__repr__()) # logb is not present in mpmath
		self.assertEqual(mne.rint(r),round(r))
		self.assertTrue((mne.round(r)==round(r)) or (r%1==0.5)) # ignore rounding 0.5 up or down.
		self.checkRelativeError(mne.tgamma(r),mpmath.gamma(r),functionName="tgamma")
		self.assertEqual(mne.trunc(r),int(r))

		self.checkRelativeError(mne.fabs(r),abs(r),functionName="fabs")

		pair = mne.frexp(abs(r))
		self.checkRelativeError(abs(r),pair[0]*mpmath.power(2,pair[1]),functionName="frexp")

		pair = mne.modf(abs(r))
		self.checkRelativeError(pair[0],(abs(r))%1,functionName="modf")
		self.assertEqual(pair[1],int(abs(r)))

		#self.assertEqual(mne.frexp(abs(r)),mne.frexp_c_test(abs(r)))
		#self.assertEqual(mne.modf(abs(r)),mne.modf_c_test(abs(r)))

		if(r==0):
			self.assertEqual(mne.sgn(r),0)
			self.assertEqual(mne.sign(r),0)
		if(r>0):
			self.assertEqual(mne.sgn(r),1)
			self.assertEqual(mne.sign(r),1)
		if(r<0):
			self.assertEqual(mne.sgn(r),-1)
			self.assertEqual(mne.sign(r),-1)

		self.checkCgalNumTraits(r)

	def checkCgalNumTraits(self,r):
		if(mne.testCgalNumTraits==False):
			print("Skipping test of CgalNumTraits")
			return
		self.assertEqual(mne.CGAL_Is_valid(r),True)
		self.checkRelativeError(mne.CGAL_Square(r),mpmath.power(r,2),functionName="pow")
		self.checkRelativeError(mne.CGAL_Sqrt(abs(r)),mpmath.sqrt(abs(r)),functionName="sqrt")
		for kk in range(5):
			k=kk+1
			self.checkRelativeError(mne.CGAL_Kth_root(k,abs(r)),mpmath.power(abs(r),1/mpmath.mpf(k)),functionName="pow")
		# CGAL uses double for intervals
		interval = mne.CGAL_To_interval(r)
		self.checkRelativeError(r,interval[0],1e-14)
		self.checkRelativeError(r,interval[1],1e-14)
		self.assertEqual(mne.CGAL_Is_finite(r),True)
		if(r==0): self.assertEqual(mne.CGAL_Sgn(r),0)
		if(r> 0): self.assertEqual(mne.CGAL_Sgn(r),1)
		if(r< 0): self.assertEqual(mne.CGAL_Sgn(r),-1)
		self.assertEqual(mne.CGAL_Sgn(0),0)
		self.assertEqual(mne.CGAL_Sgn(2.5),1)
		self.assertEqual(mne.CGAL_Sgn(-2.3),-1)

	def twoArgMathCheck(self,r1,r2):
		self.checkRelativeComplexError(mne.sin (mpmath.mpc(r1,r2)),mpmath.sin (mpmath.mpc(r1,r2)),functionName="csin")
		self.checkRelativeComplexError(mne.sinh(mpmath.mpc(r1,r2)),mpmath.sinh(mpmath.mpc(r1,r2)),functionName="csinh")
		self.checkRelativeComplexError(mne.cos (mpmath.mpc(r1,r2)),mpmath.cos (mpmath.mpc(r1,r2)),functionName="ccos")
		self.checkRelativeComplexError(mne.cosh(mpmath.mpc(r1,r2)),mpmath.cosh(mpmath.mpc(r1,r2)),functionName="ccosh")
		self.checkRelativeComplexError(mne.tan (mpmath.mpc(r1,r2)),mpmath.tan (mpmath.mpc(r1,r2)),functionName="ctan")
		self.checkRelativeComplexError(mne.tanh(mpmath.mpc(r1,r2)),mpmath.tanh(mpmath.mpc(r1,r2)),functionName="ctanh")

		self.checkRelativeComplexError(mne.exp(mpmath.mpc(r1,r2)),mpmath.exp(mpmath.mpc(r1,r2)),functionName="cexp")
		self.checkRelativeComplexError(mne.log(mpmath.mpc(r1,r2)),mpmath.log(mpmath.mpc(r1,r2)),functionName="clog")

		self.checkRelativeComplexError(mne.abs  (mpmath.mpc(r1,r2)),abs(mpmath.mpc(r1,r2)),functionName="cabs")
		self.checkRelativeComplexError(mne.conj (mpmath.mpc(r1,r2)),mpmath.conj(mpmath.mpc(r1,r2)),functionName="cconj")
		self.checkRelativeComplexError(mne.real (mpmath.mpc(r1,r2)),r1,functionName="creal")
		self.checkRelativeComplexError(mne.imag (mpmath.mpc(r1,r2)),r2,functionName="cimag")

		self.checkRelativeError(mne.atan2(r1,r2),mpmath.atan2(r1,r2),functionName="atan2")
		self.checkRelativeError(mne.fmod(abs(r1),abs(r2)),mpmath.fmod(abs(r1),abs(r2)),functionName="fmod")
		self.checkRelativeError(mne.hypot(r1,r2),mpmath.hypot(r1,r2),functionName="hypot")
		self.checkRelativeError(mne.max(r1,r2),max(r1,r2),functionName="max")
		self.checkRelativeError(mne.min(r1,r2),min(r1,r2),functionName="min")
		self.checkRelativeError(mne.pow(abs(r1),r2),mpmath.power(abs(r1),r2),functionName="pow")
		self.checkRelativeError(mne.remainder(abs(r1),abs(r2)),abs(r1)-round(abs(r1)/abs(r2))*abs(r2),functionName="remainder")
		pair = mne.remquo(abs(r1),abs(r2))
		self.checkRelativeError(pair[0],abs(r1)-round(abs(r1)/abs(r2))*abs(r2),functionName="remquo")
		self.assertEqual(pair[1]%8, round(abs(r1/r2))%8)

		self.checkRelativeError(mne.ldexp(r1,int(r2)),mpmath.mpf(r1)*mpmath.power(2,int(r2)),functionName="ldexp")

	def threeArgMathCheck(self,r1,r2,r3):
		self.checkRelativeError(mne.fma(r1,r2,r3),(mpmath.mpf(r1)*r2)+r3,functionName="fma")

	def testMathFunctions(self):
		self.assertEqual(mne.defprec , self.bits )
		zz=mpmath.acos(0)
		#print(zz.__repr__())
		#print("zz:",hex(id(zz)))
		#print("mpmath:",hex(id(mpmath)))
		a=mne.Var()
		a.val=zz
		self.assertEqual(mpmath.mp.dps , self.digs0+1 )
		#print("---- a.val=",a.val.__repr__())
		#print("---- zz   =",zz   .__repr__())
		#print("---- DPS  =",mpmath.mp.dps)
		#print("---- abs  =",abs(mpmath.mpf(a.val-zz)))
		#print("---- 10** =",self.tolerance)
		self.checkRelativeError(a.val,zz)
		self.assertEqual(mne.IsInteger, 0 )
		self.assertEqual(mne.IsSigned, 1 )
		self.assertEqual(mne.IsComplex,  0)
		if(self.bits >= 64):
			self.assertEqual(mne.RequireInitialization, 1 )
		else:
			self.assertEqual(mne.RequireInitialization, 0 )
		self.assertGreaterEqual(mne.ReadCost, 1)
		self.assertGreaterEqual(mne.AddCost, 1)
		self.assertGreaterEqual(mne.MulCost, 1)
		self.checkRelativeError(mne.highest(),self.maxval,2.1)
		self.checkRelativeError(-mne.lowest(),self.maxval,2.1)
		self.checkRelativeError(mne.Pi(),mpmath.pi)
		self.checkRelativeError(mne.Euler(),mpmath.euler)
		self.checkRelativeError(mne.Log2(),mpmath.log(2))
		self.checkRelativeError(mne.Catalan(),mpmath.catalan)
		#print("mne.epsilon() ",mne.epsilon(),"  self.expectedEpsilon = ",self.expectedEpsilon)
		self.checkRelativeError(mne.epsilon(),self.expectedEpsilon,10)
		if(self.digs0 == 6): # exception for float
			self.assertLessEqual(mne.dummy_precision(),10e-6)
		else:
			self.checkRelativeError(mpmath.log(mne.dummy_precision()/mne.epsilon())/mpmath.log(10) , mpmath.mpf(self.digs0)/10 , 1.5 )
		for x in range(50):
			if(self.nonBoostMPFR): # this looks like a bug in /usr/include/eigen3/unsupported/Eigen/MPRealSupport !
				self.assertLessEqual(abs(mne.random()-0.5),0.5)
			else:
				self.assertLessEqual(abs(mne.random()    ),1.0)
		for aa in range(4):
			for bb in range(4):
				a = (aa-3)*5
				b = bb*10
				r = mne.random(a,a+b+1)
				r2= mne.random(a,a+b+1)
				r3= mne.random(a,a+b+5)
				#print("random=",r)
				self.assertLessEqual(r,a+b+1)
				self.assertGreaterEqual(r,a)
				self.assertFalse(mne.isMuchSmallerThan(r,1,mne.epsilon()))
				self.assertTrue(mne.isMuchSmallerThan(self.expectedEpsilon,1+abs(r),mne.epsilon()))
				self.assertTrue(mne.isEqualFuzzy(r+self.expectedEpsilon*0.01,r,mne.epsilon()))
				self.checkRelativeError(mne.toLongDouble(r),float(r), 1e-14) # FIXME - should be 1e-17, but python does not support that
				self.checkRelativeError(mne.toDouble(r),float(r), 1e-14)
				self.checkRelativeError(mne.toDouble(r),float(r), 1e-14)
				self.assertEqual(mne.toLong(r),int(r))
				self.assertEqual(mne.toInt(r),int(r))
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
		mne.testArray()

	def testBasicVariable(self):
		a=mne.Var()
		self.checkRelativeError(a.val,-71.23,0.01)
		a.val=10
		self.checkRelativeError(a.val,10)
		self.checkRelativeComplexError(a.cpl,-71.23+33.23j,0.01)
		a.cpl=mpmath.mpc("1","-1")
		self.checkRelativeComplexError(a.cpl,1-1j,1e-15)
		self.checkRelativeComplexError(a.cpl,mpmath.mpc("1","-1"))

	def thisTestsExceptionReal(self):
		a=mne.Var()
		a.val="13123-123123*123"

	def thisTestsExceptionComplex(self):
		a=mne.Var()
		a.cpl="13123-123123*123-50j"

	def testWrongInput(self):
		if(self.nonBoostMPFR): # this looks like another bug in /usr/include/mpreal.h
			print("skipping this test for non-boost /usr/include/mpreal.h")
			return
		# depending on backed Real use it throws TypeError or RuntimeError
		self.assertRaises(Exception,self.thisTestsExceptionReal)
		self.assertRaises(Exception,self.thisTestsExceptionComplex)

