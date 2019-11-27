#pragma once

#include <iostream>
#include <fstream>
#include "RegularTriangulation.h"

namespace yade { // Cannot have #include directive inside.

namespace CGT {

using std::endl;
using std::cout;
using std::cerr;
using std::string;

#define NORMALIZE(vecteur) ((vecteur) = (vecteur)*(1.0/sqrt(pow((vecteur)[0],2)+pow((vecteur)[1],2)+pow((vecteur)[2],2))))

class Tens;
class Tenseur3;
class Tenseur_sym3;

CVector operator* ( Tens& tens, CVector& vect );
CVector& NormalizedCVector ( CVector& vect );


void Tenseur_produit ( CVector &v1, CVector &v2, Tenseur3 &result );
void Somme ( Tenseur3 &result, CVector &v1, CVector &v2 );

std::ostream& operator<< ( std::ostream& os,const Tenseur3& T );
std::ostream& operator<< ( std::ostream& os,const Tenseur_sym3& T );

class Tens
{
	public:
		Tens ( void ) {}
		virtual ~Tens ( void ) {}
		virtual Real operator() ( int /*i*/, int /*j*/ ) const {return 0;}
		Real Norme2 ( void );
		Real Norme ( void ) {return sqrt ( Norme2() );}
		Real Trace ( void )
		{
			return this->operator () ( 1,1 )
				   + this->operator () ( 2,2 )
				   + this->operator () ( 3,3 );
		}
};

class Tenseur3 : public Tens
{
	private:
		Real T [3] [3];

	public:
		Tenseur3 ( bool init = true );// Sp�cifier "false" pour �conomiser le temps d'initialisation du tableau
		virtual ~Tenseur3 ( void );
		Tenseur3 ( const Tenseur3& source );
		Tenseur3 ( Real a11, Real a12, Real a13,
				   Real a21, Real a22, Real a23,
				   Real a31, Real a32, Real a33 );

		Tenseur3& operator= ( const Tenseur3& source );
		Tenseur3& operator/= ( Real d );
		Tenseur3& operator+= ( const Tenseur3& source );
		Real operator() ( int i, int j ) const;
		Real &operator() ( int i, int j );

		virtual void reset ( void ) {for ( int i=0; i<3; i++ ) for ( int j=0; j<3; j++ ) T[i][j] = 0;}

};

class Tenseur_sym3 : public Tens
{
	private:
		Real T [6];

	public:
		Tenseur_sym3 ( bool init = true );// Sp�cifier "false" pour �conomiser le temps d'initialisation du tableau
		~Tenseur_sym3 ( void );
		Tenseur_sym3 ( const Tenseur_sym3& source );
		Tenseur_sym3 ( const Tenseur3& source );
		Tenseur_sym3 ( Real a11, Real a22, Real a33,
					   Real a12, Real a13, Real a23 );

		Tenseur_sym3& operator= ( const Tenseur_sym3& source );
		Tenseur_sym3& operator/= ( Real d );
		Tenseur_sym3 Deviatoric ( void ) const; //retourne la partie d�viatoire
		Real operator() ( int i, int j ) const;
		Real &operator() ( int i, int j );

		void reset ( void ) {for ( int i=0; i<6; i++ ) T[i] = 0;}

};

static const Tenseur3 NULL_TENSEUR3 ( 0,0,0,0,0,0,0,0,0 );

} // namespace CGT

} // namespace yade

