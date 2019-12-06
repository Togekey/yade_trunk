#!/bin/bash
# autopkgtest check for minieigen
# (C) 2015 Anton Gladky <gladk@debian.org>
# (C) 2019 Janek Kozicki

set -e

LIBTOTEST=${1}
DEC_DIGITS=${2}

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" 0 INT QUIT ABRT PIPE TERM
cp ${LIBTOTEST}.so ${WORKDIR}
cd ${WORKDIR}

cat <<EOF > minieigentest.py
import unittest, math, sys, mpmath
import ${LIBTOTEST} as mne

class SimpleTests(unittest.TestCase):
	def setUp(self):
		mpmath.mp.dps=${DEC_DIGITS}+1
		# tolerance = 1.001×10⁻ᵈ⁺¹, where ᵈ==${DEC_DIGITS}
		# so basically we store one more decimal digit, and expect one less decimal digit. That amounts to ignoring one (actually two) least significant digits.
		self.tolerance=(mpmath.mpf(10)**(-${DEC_DIGITS}+1))*mpmath.mpf("1.001")
	def checkRelativeError(self,a,b):
		self.assertLessEqual(abs( (mpmath.mpf(a)-mpmath.mpf(b))/mpmath.mpf(b) ),self.tolerance)
	def testVector2i(self):
		a2i = mne.Vector2i(2,1)
		b2i = mne.Vector2i(3,5)
		c2i = a2i + b2i

		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( c2i[0] , mpmath.mpf("5") )
		self.checkRelativeError( c2i[1] , mpmath.mpf("6") )

		c2i *= 2

		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( c2i[0] , mpmath.mpf("10") )
		self.checkRelativeError( c2i[1] , mpmath.mpf("12") )

	def testVector3(self):
		a3r = mne.Vector3(2.1,1.1,4.3)
		b3r = mne.Vector3(3.1,5.1,5.2)
		c3r = a3r + b3r

		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( c3r[0] , mpmath.mpf("5.2") )
		self.checkRelativeError( c3r[1] , mpmath.mpf("6.2") )
		self.checkRelativeError( c3r[2] , mpmath.mpf("9.5") )

		c3r *= 3

		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( c3r[0] , mpmath.mpf("15.6") )
		self.checkRelativeError( c3r[1] , mpmath.mpf("18.6") )
		self.checkRelativeError( c3r[2] , mpmath.mpf("28.5") )

	def testMatrix3Test(self):
		ll=1
		print("------------------------------------------------------",ll);ll+=1;
		a3m=mne.Matrix3(1,2,3,
		                4,5,6,
		                7,8,9)
		print("------------------------------------------------------",ll);ll+=1;
		b3m=a3m.transpose()
		print("------------------------------------------------------",ll);ll+=1;
		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		print("------------------------------------------------------",ll);ll+=1;
		aaaaa=b3m[0][0]
		print("------------------------------------------------------",ll,aaaaa);ll+=1;
		print(b3m[0][0])
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( b3m[0][0] , mpmath.mpf("1") )
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( b3m[0][1] , mpmath.mpf("4") )
		self.checkRelativeError( b3m[0][2] , mpmath.mpf("7") )
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( b3m[1][0] , mpmath.mpf("2") )
		self.checkRelativeError( b3m[1][1] , mpmath.mpf("5") )
		self.checkRelativeError( b3m[1][2] , mpmath.mpf("8") )
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( b3m[2][0] , mpmath.mpf("3") )
		self.checkRelativeError( b3m[2][1] , mpmath.mpf("6") )
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( b3m[2][2] , mpmath.mpf("9") )

		c3m=a3m.diagonal()
		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( c3m[0] , mpmath.mpf("1") )
		self.checkRelativeError( c3m[1] , mpmath.mpf("5") )
		self.checkRelativeError( c3m[2] , mpmath.mpf("9") )

		self.checkRelativeError( a3m.maxAbsCoeff() , mpmath.mpf("9") )


	def testQuaternion(self):
		ll=1
		print("------------------------------------------------------",ll);ll+=1;
		q1 = mne.Quaternion.Identity
		print("------------------------------------------------------",ll,q1);ll+=1;
		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( q1[3] , mpmath.mpf("1") )
		print("------------------------------------------------------",ll);ll+=1;

		q2 = q1.inverse()
		self.checkRelativeError( q2[3] , mpmath.mpf("1") )
		print("------------------------------------------------------",ll);ll+=1;
		q3=mne.Quaternion(axis=mne.Vector3(1,0,0),angle=mpmath.pi/2.0)
		m3q=q3.toRotationMatrix()
		self.checkRelativeError( m3q[0][0] , mpmath.mpf("1") )
		print(m3q[1][2].__repr__())
		print("------------------------------------------------------",ll);ll+=1;
		self.checkRelativeError( m3q[1][2] , mpmath.mpf("-1") )
		print(q3)
		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )

		print("------------------------------------------------------",ll);ll+=1;
		q4 = mne.Quaternion.Identity
		q4.setFromTwoVectors(mne.Vector3(1,2,3),mne.Vector3(2,3,4))
		print(q4.norm().__repr__())
		print("------------------------------------------------------",ll);ll+=1;
		self.assertEqual(mpmath.mp.dps , ${DEC_DIGITS}+1 )
		self.checkRelativeError( q4.norm() , mpmath.mpf("1") )
		print("------------------------------------------------------",ll);ll+=1;

if __name__ == '__main__':
		unittest.main(testRunner=unittest.TextTestRunner(stream=sys.stdout, verbosity=2))
EOF

function handle_error() {
	ls -la
	gdb --batch -ex "bt full" `which python3` ./core
	exit 1
}

echo 'Test Python3'
ulimit -c unlimited
python3 minieigentest.py || handle_error
echo "Python3 run: OK"

cd
rm -rf ${WORKDIR}

