# -*- coding: utf-8 -*-
# an extended version of autopkgtest check for minieigen
# (C) 2015 Anton Gladky <gladk@debian.org>
# (C) 2019 Janek Kozicki

import unittest, math, sys
import minieigen as mne
import yade

import testMathHelper as mpmath
from   testMathHelper import mpc

class ExtendedMinieigenTests(unittest.TestCase):
	def setUp(self):
		self.digs1=yade.config.highPrecisionDecimalPlaces+1
		mpmath.mp.dps=self.digs1
		# tolerance = 1.001×10⁻ᵈ⁺¹, where ᵈ==yade.config.highPrecisionDecimalPlaces
		# so basically we store one more decimal digit, and expect one less decimal digit. That amounts to ignoring one (actually two) least significant digits.
		self.tolerance=(mpmath.mpf(10)**(-yade.config.highPrecisionDecimalPlaces+1))*mpmath.mpf("1.001")
	def checkRelativeError(self,a,b):
		if b!=0:
			self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/mpmath.mpf(b) ),self.tolerance)
		else:
			self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/self.tolerance ),self.tolerance)
	def checkRelativeComplexError(self,a,b):
		self.assertLessEqual(abs( (mpmath.mpc(a)-mpmath.mpc(b))/mpmath.mpc(b) ),self.tolerance)

	def testMpmath(self):
		self.assertEqual(2    ,float(mpmath.mpf(2)))
		self.assertEqual(2    ,float(mpmath.mpf("2")))
		self.assertEqual(2**3 ,mpmath.mpf("2")**3)
		self.assertEqual(2/4  ,mpmath.mpf("2")/4)
		self.assertEqual(3/2,3/mpmath.mpf("2"))
		self.assertEqual(3*2,3*mpmath.mpf("2"))
		self.assertEqual(2-3  ,mpmath.mpf("2")-3)

		self.assertEqual(2    ,complex(mpmath.mpc(2)))
		self.assertEqual(2+5j ,complex(mpmath.mpc(2+5j)))
		self.assertEqual(2    ,complex(mpmath.mpc("2")))
		self.assertEqual(2+5j ,complex(mpmath.mpc("2","5")))
		self.assertEqual(2-3  ,mpmath.mpc("2")-3)
		self.assertEqual(2/4  ,mpmath.mpc("2")/4)
		self.assertEqual(3/2,3/mpmath.mpc("2"))
		self.assertEqual(5,abs(mpmath.mpc("-3","-4")))


	def testVector2i(self):
		a2i = mne.Vector2i(2,1)
		b2i = mne.Vector2i(3,5)
		c2i = a2i + b2i
		#self.assertEqual(c2i.eigenFlags(),352)
		#self.assertEqual(c2i.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c2i[0] , mpmath.mpf("5") )
		self.checkRelativeError( c2i[1] , mpmath.mpf("6") )

		c2i *= 2

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c2i[0] , mpmath.mpf("10") )
		self.checkRelativeError( c2i[1] , mpmath.mpf("12") )

		from minieigen import Vector2i
		self.assertEqual( c2i , eval(c2i.__str__()) )

	def testVector2(self):
		a2 = mne.Vector2(2,1)
		b2 = mne.Vector2(3,5)
		c2 = a2 + b2
		#self.assertEqual(c2.eigenFlags(),352)
		#self.assertEqual(c2.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c2[0] , mpmath.mpf("5") )
		self.checkRelativeError( c2[1] , mpmath.mpf("6") )

		c2 *= 2

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c2[0] , mpmath.mpf("10") )
		self.checkRelativeError( c2[1] , mpmath.mpf("12") )

		from minieigen import Vector2
		self.assertEqual( c2 , eval(c2.__str__()) )

	def testVector2c(self):
		a2c = mne.Vector2c(2-10j,1)
		b2c = mne.Vector2c(3,5)
		c2c = a2c + b2c
		#self.assertEqual(c2c.eigenFlags(),352)
		#self.assertEqual(c2c.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( c2c[0] , mpmath.mpc("5","-10") )
		self.checkRelativeComplexError( c2c[1] , mpmath.mpf("6") )

		c2c *= 2

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( c2c[0] , mpmath.mpc("10","-20") )
		self.checkRelativeComplexError( c2c[1] , mpmath.mpf("12") )

		from minieigen import Vector2c
		self.assertEqual( c2c , eval(c2c.__str__()) )

	def testVector3i(self):
		a3i = mne.Vector3i(2,1,4)
		b3i = mne.Vector3i(3,5,5)
		c3i = a3i + b3i
		#self.assertEqual(c3i.eigenFlags(),352)
		#self.assertEqual(c3i.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3i[0] , mpmath.mpf("5") )
		self.checkRelativeError( c3i[1] , mpmath.mpf("6") )
		self.checkRelativeError( c3i[2] , mpmath.mpf("9") )

		c3i *= 3

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3i[0] , mpmath.mpf("15") )
		self.checkRelativeError( c3i[1] , mpmath.mpf("18") )
		self.checkRelativeError( c3i[2] , mpmath.mpf("27") )

		from minieigen import Vector3i
		self.assertEqual( c3i , eval(c3i.__str__()) )

	def testVector3(self):
		a3r = mne.Vector3(2.1,1.1,4.3)
		b3r = mne.Vector3(3.1,5.1,5.2)
		c3r = a3r + b3r
		#self.assertEqual(c3r.eigenFlags(),352)
		#self.assertEqual(c3r.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3r[0] , mpmath.mpf("5.2") )
		self.checkRelativeError( c3r[1] , mpmath.mpf("6.2") )
		self.checkRelativeError( c3r[2] , mpmath.mpf("9.5") )

		c3r *= 3

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3r[0] , mpmath.mpf("15.6") )
		self.checkRelativeError( c3r[1] , mpmath.mpf("18.6") )
		self.checkRelativeError( c3r[2] , mpmath.mpf("28.5") )

		from minieigen import Vector3
		self.checkRelativeError( c3r[0] , eval(c3r.__str__())[0] )
		self.checkRelativeError( c3r[1] , eval(c3r.__str__())[1] )
		self.checkRelativeError( c3r[2] , eval(c3r.__str__())[2] )

	def testVector3c(self):
		a3c = mne.Vector3c(2.1+1j,1.1+2.5j,4.3-1j)
		b3c = mne.Vector3c(3.1,5.1,5.2)
		c3c = a3c + b3c
		#self.assertEqual(c3c.eigenFlags(),352)
		#self.assertEqual(c3c.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( c3c[0] , mpmath.mpc("5.2","1") )
		self.checkRelativeComplexError( c3c[1] , mpmath.mpc("6.2","2.5") )
		self.checkRelativeComplexError( c3c[2] , mpmath.mpc("9.5","-1") )

		c3c *= 3

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( c3c[0] , mpmath.mpc("15.6","3") )
		self.checkRelativeComplexError( c3c[1] , mpmath.mpc("18.6","7.5") )
		self.checkRelativeComplexError( c3c[2] , mpmath.mpc("28.5","-3") )

		from minieigen import Vector3c
		self.checkRelativeComplexError( c3c[0] , eval(c3c.__str__())[0] )
		self.checkRelativeComplexError( c3c[1] , eval(c3c.__str__())[1] )
		self.checkRelativeComplexError( c3c[2] , eval(c3c.__str__())[2] )

	def testVector3na(self):
		if(mne.vectorize == False): return
		a3a = mne.Vector3na(2.1,1.1,4.3)
		b3a = mne.Vector3na(3.1,5.1,5.2)
		c3a = a3a + b3a
		#self.assertEqual(c3a.eigenFlags(),352)
		#self.assertEqual(c3a.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3a[0] , mpmath.mpf("5.2") )
		self.checkRelativeError( c3a[1] , mpmath.mpf("6.2") )
		self.checkRelativeError( c3a[2] , mpmath.mpf("9.5") )

		c3a *= 3

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3a[0] , mpmath.mpf("15.6") )
		self.checkRelativeError( c3a[1] , mpmath.mpf("18.6") )
		self.checkRelativeError( c3a[2] , mpmath.mpf("28.5") )

		from minieigen import Vector3na
		self.checkRelativeError( c3a[0] , eval(c3a.__str__())[0] )
		self.checkRelativeError( c3a[1] , eval(c3a.__str__())[1] )
		self.checkRelativeError( c3a[2] , eval(c3a.__str__())[2] )

	def testVector4(self):
		# The Vector4 bug was fixed only recently, don't test if there's nothing to test
		if(not hasattr(mne, 'Vector4')): return
		a4r = mne.Vector4(2.1,1.1,4.3,5.5)
		b4r = mne.Vector4(3.1,5.1,5.2,-5.0)
		c4r = a4r + b4r
		#self.assertEqual(c4r.eigenFlags(),352)
		#self.assertEqual(c4r.eigenStorageOrder(),0)

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c4r[0] , mpmath.mpf("5.2") )
		self.checkRelativeError( c4r[1] , mpmath.mpf("6.2") )
		self.checkRelativeError( c4r[2] , mpmath.mpf("9.5") )
		self.checkRelativeError( c4r[3] , mpmath.mpf("0.5") )

		c4r *= 3

		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c4r[0] , mpmath.mpf("15.6") )
		self.checkRelativeError( c4r[1] , mpmath.mpf("18.6") )
		self.checkRelativeError( c4r[2] , mpmath.mpf("28.5") )
		self.checkRelativeError( c4r[3] , mpmath.mpf("1.5") )

		from minieigen import Vector4
		self.checkRelativeError( c4r[0] , eval(c4r.__str__())[0] )
		self.checkRelativeError( c4r[1] , eval(c4r.__str__())[1] )
		self.checkRelativeError( c4r[2] , eval(c4r.__str__())[2] )
		self.checkRelativeError( c4r[3] , eval(c4r.__str__())[3] )

	def testMatrix3Test(self):
		a3m=mne.Matrix3(1,2,3,
		                4,5,6,
		                7,8,9)
		#self.assertEqual(a3m.eigenFlags(),352)
		#self.assertEqual(a3m.eigenStorageOrder(),0)
		b3m=a3m.transpose()
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( b3m[0][0] , mpmath.mpf("1") )
		self.checkRelativeError( b3m[0][1] , mpmath.mpf("4") )
		self.checkRelativeError( b3m[0][2] , mpmath.mpf("7") )
		self.checkRelativeError( b3m[1][0] , mpmath.mpf("2") )
		self.checkRelativeError( b3m[1][1] , mpmath.mpf("5") )
		self.checkRelativeError( b3m[1][2] , mpmath.mpf("8") )
		self.checkRelativeError( b3m[2][0] , mpmath.mpf("3") )
		self.checkRelativeError( b3m[2][1] , mpmath.mpf("6") )
		self.checkRelativeError( b3m[2][2] , mpmath.mpf("9") )

		c3m=a3m.diagonal()
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( c3m[0] , mpmath.mpf("1") )
		self.checkRelativeError( c3m[1] , mpmath.mpf("5") )
		self.checkRelativeError( c3m[2] , mpmath.mpf("9") )

		self.checkRelativeError( a3m.maxAbsCoeff() , mpmath.mpf("9") )

		from minieigen import Matrix3
		for i in range(3):
			for j in range(3):
				self.checkRelativeError( b3m[i][j] , eval(b3m.__str__())[i][j] )
		#print(b3m.__str__())

	def testMatrix3cTest(self):
		a3m=mne.Matrix3c(1+1j,2,3,
		                4,5,6,
		                7,8,9-9j)
		#self.assertEqual(a3m.eigenFlags(),352)
		#self.assertEqual(a3m.eigenStorageOrder(),0)
		b3m=a3m.transpose()
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( b3m[0][0] , mpmath.mpc("1","1") )
		self.checkRelativeComplexError( b3m[0][1] , mpmath.mpc("4","0") )
		self.checkRelativeComplexError( b3m[0][2] , mpmath.mpc("7","0") )
		self.checkRelativeComplexError( b3m[1][0] , mpmath.mpc("2","0") )
		self.checkRelativeComplexError( b3m[1][1] , mpmath.mpc("5","0") )
		self.checkRelativeComplexError( b3m[1][2] , mpmath.mpc("8","0") )
		self.checkRelativeComplexError( b3m[2][0] , mpmath.mpc("3","0") )
		self.checkRelativeComplexError( b3m[2][1] , mpmath.mpc("6","0") )
		self.checkRelativeComplexError( b3m[2][2] , mpmath.mpc("9","-9") )

		c3m=a3m.diagonal()
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeComplexError( c3m[0] , mpmath.mpc("1","1") )
		self.checkRelativeComplexError( c3m[1] , mpmath.mpc("5","0") )
		self.checkRelativeComplexError( c3m[2] , mpmath.mpc("9","-9") )

		self.checkRelativeComplexError( a3m.maxAbsCoeff() , abs(mpmath.mpc("9","-9")) )

		from minieigen import Matrix3c
		for i in range(3):
			for j in range(3):
				self.checkRelativeComplexError( b3m[i][j] , eval(b3m.__str__())[i][j] )
		#print(b3m.__str__())


	def testQuaternion(self):
		q1 = mne.Quaternion.Identity
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( q1[3] , mpmath.mpf("1") )
		#self.assertEqual(q1.eigenFlags(),32)
		#self.assertEqual(q1.eigenStorageOrder(),0)

		q2 = q1.inverse()
		self.checkRelativeError( q2[3] , mpmath.mpf("1") )
		if(yade.config.highPrecisionMpmath):
			q3=mne.Quaternion(axis=mne.Vector3(1,0,0),angle=mpmath.pi/2.0)
		else:
			q3=mne.Quaternion(axis=mne.Vector3(1,0,0),angle=1.570796326794896619231321691639)
		m3q=q3.toRotationMatrix()
		self.checkRelativeError( m3q[0][0] , mpmath.mpf("1") )
		#print(m3q[1][2].__repr__())
		self.checkRelativeError( m3q[1][2] , mpmath.mpf("-1") )
		#print(q3)
		self.assertEqual(mpmath.mp.dps , self.digs1 )

		q4 = mne.Quaternion.Identity
		q4.setFromTwoVectors(mne.Vector3(1,2,3),mne.Vector3(2,3,4))
		#print(q4.norm().__repr__())
		self.assertEqual(mpmath.mp.dps , self.digs1 )
		self.checkRelativeError( q4.norm() , mpmath.mpf("1") )

		from minieigen import Quaternion
		self.checkRelativeError( q3[0] , eval(q3.__str__())[0] )
		self.checkRelativeError( q3[1] , eval(q3.__str__())[1] )
		self.checkRelativeError( q3[2] , eval(q3.__str__())[2] )
		self.checkRelativeError( q3[3] , eval(q3.__str__())[3] )
		#print(q3.__str__())

