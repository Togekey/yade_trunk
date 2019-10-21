/*CWBoon 2015 */
#ifdef YADE_POTENTIAL_BLOCKS
#pragma once

#include <core/Shape.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <vector>
//#include <lib/base/openmp-accu.hpp>

namespace yade { // Cannot have #include directive inside.


class PotentialBlock : public Shape {
public:
	struct Planes {
		vector<int> vertexID;
	};
	//		struct Vertices{ vector<int> edgeID;  vector<int> planeID; };
	//		struct Edges{ vector<int> vertexID; };

	void addPlaneStruct();
	//		void addVertexStruct();
	//		void addEdgeStruct();

	vector<Planes> planeStruct;
	//		vector<Vertices> vertexStruct;
	//		vector<Edges> edgeStruct;

	//		Eigen::MatrixXd Amatrix;
	//		Eigen::MatrixXd Dmatrix;
	virtual ~PotentialBlock();
	void postLoad(PotentialBlock&);

	double getDet(const Eigen::MatrixXd A);
	double getSignedArea(const Vector3r pt1, const Vector3r pt2, const Vector3r pt3);
	void   calculateVertices();
	void   calculateInertia(Vector3r& centroid, Real& Ixx, Real& Iyy, Real& Izz, Real& Ixy, Real& Ixz, Real& Iyz);

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(PotentialBlock,Shape,"Geometry of PotentialBlock.",
		((bool, isLining, false,, "Whether particle is part of tunnel lining (used in the RockLining.cpp script)"))
		((double, liningStiffness, pow(10.0,8),, "Lining stiffness"))
		((double, liningFriction, 20.0,, "Lining friction"))
		((double, liningLength, 0.0,, "Lining spacing between nodes"))
		((double, liningTensionGap, 0.0,, "Numerical gap between lining and block to allowing tension to be calculated"))
		((Vector3r, liningNormalPressure, Vector3r(0,0,0),, "Normal pressure acting on lining"))
		((Vector3r, liningTotalPressure, Vector3r(0,0,0),, "Total pressure acting on lining"))
		((bool, isBoundary, false,, "Whether the particle is part of a boundary block"))
//			((bool, isEastBoundary, false,Attr::hidden, "Boolean used for a case history"))
		((bool, isBolt, false,, "Whether a block is part of a bolt (used in the Rockbolt.cpp script)"))
		((bool, fixedNormal, false,, "Whether to fix the contact normal at a boundary, using boundaryNormal"))
		((Vector3r, boundaryNormal, Vector3r::Zero(),,"Normal direction of boundary if fixedNormal=True"))
			((bool, AabbMinMax, false,, "Whether the exact Aabb should be calculated. If false, an approximate cubic Aabb is defined with edges of ``2R``"))
		((Vector3r, minAabb, Vector3r::Zero(),,"Min from box centre: Used for visualisation in vtk"))
		((Vector3r, maxAabb, Vector3r::Zero(),,"Max from box centre: Used for visualisation in vtk"))
//			((Vector3r, minAabbRotated, Vector3r::Zero(),,"(not used: Actually used in the BlockGen, but not needed!) min from box centre"))
//			((Vector3r, maxAabbRotated, Vector3r::Zero(),,"(not used: Actually used in the BlockGen, but not needed!) max from box centre"))
//			((Vector3r, halfSize, Vector3r::Zero(),Attr::hidden,"max from box centre"))
//			((vector<Vector3r>, node, ,Attr::hidden, "(not used)"))
//			((Real, gridVol, ,Attr::hidden, "(not used)"))
//			((Quaternionr, oriAabb, Quaternionr::Identity(),, "(not used) original aabb "))
		((Real, r, 0.0,, "r in Potential Particles"))
		((Real, R, 0.0,, "R in Potential Particles. If left zero, a default value is calculated as half the distance of the farthest vertices"))
		((Real, k, 0.0,, "k in Potential Particles (not used)"))
		((Real, volume, ,, "Volume"))
		((Vector3r, inertia, Vector3r::Zero(),, "Principal inertia tensor"))
		((Quaternionr, orientation, Quaternionr::Identity(),, "Principal orientation"))
		((int, id, -1,, "Particle id (for graphics in vtk output)")) //TODO: Check if we can use the body id instead in all instances and delete this attribute
		((bool, erase, false,, "Parameter to mark particles to be removed (for excavation)"))
		((vector<bool>, intactRock, false,, "Property for plane")) //TODO: Provide more info on the functionality of this attribute
//			((vector<bool>, isBoundaryPlane, ,Attr::hidden, "Property for plane (not used)"))
//			((vector<double>, hwater, ,Attr::hidden, "Property for plane (not used), height of pore water"))
//			((vector<double>, JRC, ,Attr::hidden, "Property for plane (not used), rock joint"))
//			((vector<double>, JCS, ,Attr::hidden, "Property for plane (not used), rock joint"))
//			((vector<double>, asperity, ,Attr::hidden, "Property for plane (not used), rock joint"))
//			((vector<double>, sigmaC, ,Attr::hidden,"Property for plane (not used), rock joint"))
		((vector<double>, phi_b, ,, "Basic friction angle of each face (property for plane, rock joint)"))
		((vector<double>, phi_r, ,, "Residual friction angle of each face (property for plane, rock joint)"))
		((vector<double>, cohesion, ,, "Cohesion (stress) of each face (property for plane, rock joint)"))
		((vector<double>, tension, ,, "Tension (stress) of each face (property for plane, rock joint)"))
		((vector<int>, jointType, ,, "jointType"))
//			((vector<double>, lambda0, ,Attr::hidden, "Property for plane (not used), heat"))
//			((vector<double>, kn, ,Attr::hidden, "Property for plane, rock joint (not used: would be used in that each face could have it's own stiffness properties (not developed))"))
//			((vector<double>, ks, ,Attr::hidden, "Property for plane, rock joint (not used: would be used in that each face could have it's own stiffness properties (not developed))"))
//			((vector<double>, heatCapacity, ,Attr::hidden, "Property for plane, rock joint"))
//			((vector<double>, rFactor, ,Attr::hidden, "(not used), individual factor for r"))
		((vector<Vector3r>, vertices,,,"Vertices"))
		//((Eigen::MatrixXd , Amatrix, ,, "a "))
		//((Eigen::MatrixXd , Dmatrix, ,, "b "))
//			((double, waterVolume, ,, "volume of body submerged in water"))
//			((vector<Vector3r>, verticesCD, ,, "vertices"))
		((vector<double>, a, ,, "List of a coefficients of plane normals"))
		((vector<double>, b, ,, "List of b coefficients of plane normals"))
		((vector<double>, c, ,, "List of c coefficients of plane normals"))
		((vector<double>, d, ,, "List of d coefficients of plane equations"))
		, /*ctor*/
		createIndex();
		#if 0
		for (int i=0; i<a.size(); i++){
			Amatrix(i,0) = a[i]; Amatrix(i,1)=b[i]; Amatrix(i,2)=c[i];
			Dmatrix(i,0) = d[i] + r;
	 	}
		#endif
	);
	// clang-format on
	//#endif
	REGISTER_CLASS_INDEX(PotentialBlock, Shape);
};

REGISTER_SERIALIZABLE(PotentialBlock);

} // namespace yade

#ifdef __cplusplus
extern "C" {
#endif
void dgesv_(const int* N, const int* nrhs, double* Hessian, const int* lda, int* ipiv, double* gradient, const int* ldb, int* info);
void dsyev_(const char* jobz, const char* uplo, const int* N, double* A, const int* lda, double* W, double* work, int* lwork, int* info);
#ifdef __cplusplus
};
#endif
#endif // YADE_POTENTIAL_BLOCKS
