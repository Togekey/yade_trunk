/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              * 
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//#ifndef FLOW_GUARD
//#define FLOW_GUARD

#ifdef PARTIALSAT
#include "FlowEngine_PartialSatClayEngineT.hpp"
#include<Eigen/SparseLU>
#include<core/PartialEngine.hpp>
#include<core/State.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<pkg/common/Dispatching.hpp>
#include<core/Scene.hpp>
#include<core/Omega.hpp>

#ifdef FLOW_ENGINE
//#include<pkg/pfv/FlowEngine.hpp>
#include<lib/triangulation/Tesselation.h>
#include<lib/triangulation/FlowBoundingSphere.hpp>
#include "FlowEngine_FlowEngineT.hpp"
#include<pkg/dem/TesselationWrapper.hpp>
#include<lib/triangulation/Network.hpp>
#endif

#ifdef CHOLMOD_LIBS
	#include <cholmod.h>
#endif

class PartialSatCellInfo : public FlowCellInfo_PartialSatClayEngineT
{
	public:
	double saturation;//the saturation of single pore (will be used in quasi-static imbibition and dynamic flow)
	double porosity;
	//double solidLine [4][4];//the length of intersecting line between sphere and facet. [i][j] is for facet "i" and sphere (facetVertices)"[i][j]". Last component [i][3] for 1/sumLines in the facet "i" (used by chao).
	double dsdp; // the change of saturation for given capillary pressure 
	bool crack;
	Real crackArea;
	bool isExposed;
	//DynamicTwoPhaseFlow 
	//std::vector<double> entryPressure;
	//std::vector<double> entrySaturation;

	PartialSatCellInfo (void)
	{
		saturation = 1.0;
		porosity=1.0;
		dsdp = 0;
		crack=false;
		isExposed=false; // flag for determining if a pore is exposed to atmosphere, which controls the pressure force calculations
	}

	inline double& sat (void) {return saturation;}
	
};

class PartialSatVertexInfo : public FlowVertexInfo_PartialSatClayEngineT {
	public:
	//same here if needed
};

typedef CGT::_Tesselation<CGT::TriangulationTypes<PartialSatVertexInfo,PartialSatCellInfo> > PartialSatTesselation;
#ifdef LINSOLV
#define PartialSatBoundingSphere CGT::PartialSatLinSolv<PartialSatTesselation>
//class PartialSatBoundingSphere; // : public CGT::FlowBoundingSphereLinSolv<PartialSatTesselation> {};
#endif

typedef TemplateFlowEngine_PartialSatClayEngineT<PartialSatCellInfo,PartialSatVertexInfo,PartialSatTesselation,PartialSatBoundingSphere> PartialSatClayEngineT;

REGISTER_SERIALIZABLE(PartialSatClayEngineT);
YADE_PLUGIN((PartialSatClayEngineT));
class PartialSatClayEngine : public PartialSatClayEngineT
{

	public:			
		//typedef TemplateFlowEngine_FlowEngineT<FlowCellInfo_FlowEngineT,FlowVertexInfo_FlowEngineT> FlowEngineT;			
		typedef PartialSatClayEngineT::Tesselation					Tesselation;
		typedef PartialSatClayEngineT::RTriangulation					RTriangulation;
		typedef PartialSatClayEngineT::FiniteCellsIterator				FiniteCellsIterator;
		typedef PartialSatClayEngineT::CellHandle						CellHandle;
		typedef PartialSatClayEngineT::VertexHandle	VertexHandle;
		typedef std::vector<CellHandle>		VectorCell;
		typedef typename VectorCell::iterator		VCellIterator;
	public :
	double dsdp(CellHandle& cell);
	void initializeSaturations(FlowSolver& flow);
	void setSaturationFromPcS(CellHandle& cell);
	void setCellsDSDP(FlowSolver& flow);
	void updateSaturation(FlowSolver& flow);
//	void triangulate(FlowSolver& flow);
	double diagonalSaturationContribution(CellHandle cell);
	double RHSSaturationContribution(CellHandle cell);	
	void collectParticleSuction(FlowSolver& flow);
	void swellParticles();
	void setOriginalParticleValues();
	void exposureRecursion(CellHandle cell);
	void determineFracturePaths();
	void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool big, bool dynamic );
	void insertMicroPores(const float fracMicroPore);
	bool findInscribedRadiusAndLocation(CellHandle& cell, std::vector<double>& coordAndRad);
	bool checkSphereContainedInTet(CellHandle& cell, std::vector<double>& coordAndRad);
	// fracture network 
	void trickPermeability(Solver* flow);
	void interpolateCrack(Tesselation& Tes,Tesselation& NewTes);
	void computeFracturePerm(RTriangulation::Facet_circulator& facet,Real aperture,RTriangulation::Finite_edges_iterator& edge);
	void circulateFacets(RTriangulation::Finite_edges_iterator& edge,Real aperture);
//	void setPositionsBuffer(bool current);
	Real leakOffRate;
    	Real averageAperture;
	Real averageFracturePermeability;
    	Real maxAperture;	
	Real crackArea;
	Real totalFractureArea;
	

	
	virtual void initializeVolumes(FlowSolver& flow);
	virtual void buildTriangulation(double pZero, Solver& flow);
	virtual void initSolver(FlowSolver& flow);
	virtual void action();

	virtual ~PartialSatClayEngine();


	double getCellSaturation(Vector3r pos){return solver->getCellSaturation(pos[0], pos[1], pos[2]);}
	CELL_SCALAR_GETTER(Real,.sat(),cellSaturation);
	CELL_SCALAR_SETTER(double,.sat(),setCellSaturation);

	//FlowEngineT* flow;
	
//	PartialSatClayEngineT* flow;

	//We can overload every functions of the base engine to make it behave differently
	//if we overload action() like this, this engine is doing nothing in a standard timestep, it can still have useful functions
//	virtual void action() {};

	void saveUnsatVtk(const char* folder, bool withBoundaries);
//	void computeOnePhaseFlow() {scene = Omega::instance().getScene().get(); if (!solver) cerr<<"no solver!"<<endl; solver->gaussSeidel(scene->dt);initSolver(*solver);}


//	CELL_SCALAR_GETTER(bool,.isWRes,cellIsWRes)
//	CELL_SCALAR_SETTER(Real,.dvTPF,setCellDV) //Temporary function to allow for simulations in Python
	
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PartialSatClayEngine,PartialSatClayEngineT,"documentation here",
	((double,lmbda,0.2,,"Lambda parameter for Van Genuchten model"))
	((double, pAir,0,,"Air pressure for calculation of capillary pressure (Pair - Pwater)"))
	((double, Po,1.5,,"Po parameter for Van Genuchten model"))
	((double, nUnsatPerm,5,,"n parameter for empirical relative saturation based permeability relationship"))
	((double, SrM,0.01,,"residual saturation for empirical relative saturation based permeability relationship"))
	((double, SsM,1.,,"saturated saturation for empirical relative saturation based permeability relationship"))
	((bool, partialSatEngine,1,,"Activates the partial sat clay engine"))
	((bool, allCellsFractured,0,,"use to simulate all pores fractured for debugging purposes only"))
	((bool, freezeRadii,0,,"use to stop swelling debugging purposes only"))
	((bool, crackModelActive,0,,"Activates the parallel plate approximation model for facets connected to cohesionBroken edges"))
	((double,alpham,2.6048e-08,,"alpha parameter for particle volumetric strain model MPa^-1"))
	((double,betam,2.10206e-08,,"beta parameter for particle volumetric strain model MPa^-1"))
	((bool,particleSwelling,1,,"set false to neglect particle swelling"))
	((bool,freeSwelling,1,,"if true, boundary forces are computed with pAir pressure only"))
	((double,totalVolChange,0,,"tracks the total volumetric strain that occured in each step"))
	((double,matricSuctionRatio,1,,"The ratio of matric:osmotic suction. Facet forces computed for matricSuction fraction only."))
	((double,residualAperture,0.,,"residual aperture of induced cracks"))
	((double,apertureFactor,1.,,"factor to consider tortuosity"))
	((bool,calcCrackArea,true,,"The amount of crack per pore is updated if calcCrackArea=True")) 
	((double,microStructureE,1e6,,"The amount of crack per pore is updated if calcCrackArea=True")) 
	((double,microStructureNu,0.3,,"The amount of crack per pore is updated if calcCrackArea=True")) 
	((double,microStructurePhi,18.,,"The amount of crack per pore is updated if calcCrackArea=True")) 
	((double,microStructureRho,2600,,"The amount of crack per pore is updated if calcCrackArea=True")) 
	((double,microStructureAdh,6e6,,"Adhesion between microstructure particles"))
	((bool,swelling,true,,"turn just particle swelling off (for debug)"))
	((bool,suction,true,,"turn just particle suction off (for debug)"))
	((bool,volumes,true,,"turn just particle volumes off (for debug)"))
	((double,minMicroRadFrac,0.1,,"Used during sphere insertion checks, if inscribed sphere contacts facet it cannot be reduced further than minMicroRadFrac*originalInscribedRadius"))
//	((bool,saturation,true,,"turn just particle saturation off (for debug)"))
	,/*PartialSatClayEngineT()*/,
	solver = shared_ptr<FlowSolver> (new FlowSolver);
	,
	.def("setCellSaturation",&PartialSatClayEngine::setCellSaturation,(boost::python::arg("id"),boost::python::arg("temperature")),"set temperature in cell 'id'.")
	.def("getCellSaturation",&PartialSatClayEngine::getCellSaturation,(boost::python::arg("pos")),"Measure cell saturation in position pos[0],pos[1],pos[2]")
//	.def("getCellSaturation",&TwoPhaseFlowEngine::getCellSaturation,"Get saturation of cell")
	.def("saveUnsatVtk",&PartialSatClayEngine::saveUnsatVtk,(boost::python::arg("folder")="./VTK",boost::python::arg("withBoundaries")=false),"Save pressure and saturation field in vtk format. Specify a folder name for output. The cells adjacent to the bounding spheres are generated conditionally based on :yref:`withBoundaries` (not compatible with periodic boundaries)")
	.def("insertMicroPores",&PartialSatClayEngine::insertMicroPores,(boost::python::arg("fracMicroPores")),"run to inscribe spheres in a desired fraction of existing pores.")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PartialSatClayEngine);

#endif //TwoPhaseFLOW
 
