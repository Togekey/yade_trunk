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
	//DynamicTwoPhaseFlow 
	//std::vector<double> entryPressure;
	//std::vector<double> entrySaturation;

	PartialSatCellInfo (void)
	{
		saturation = 1.0;
		porosity=1.0;
		dsdp = 0;
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
	void triangulate(FlowSolver& flow);
	double diagonalSaturationContribution(CellHandle cell);
	double RHSSaturationContribution(CellHandle cell);	
	void collectParticleSaturation(FlowSolver& flow);
	void swellParticles();
	void setOriginalParticleVolumes();


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
	((bool, partialSatEngine,1,,"Activates the partial sat clay engine"))
	((double,alpham,2.6048e-08,,"alpha parameter for particle volumetric strain model MPa^-1"))
	((double,betam,2.10206e-08,,"beta parameter for particle volumetric strain model MPa^-1"))
	((bool,particleSwelling,1,,"set false to neglect particle swelling"))
	,/*PartialSatClayEngineT()*/,
	solver = shared_ptr<FlowSolver> (new FlowSolver);
	,
	.def("setCellSaturation",&PartialSatClayEngine::setCellSaturation,(boost::python::arg("id"),boost::python::arg("temperature")),"set temperature in cell 'id'.")
	.def("getCellSaturation",&PartialSatClayEngine::getCellSaturation,(boost::python::arg("pos")),"Measure cell saturation in position pos[0],pos[1],pos[2]")
//	.def("getCellSaturation",&TwoPhaseFlowEngine::getCellSaturation,"Get saturation of cell")
	.def("saveUnsatVtk",&PartialSatClayEngine::saveUnsatVtk,(boost::python::arg("folder")="./VTK",boost::python::arg("withBoundaries")=false),"Save pressure and saturation field in vtk format. Specify a folder name for output. The cells adjacent to the bounding spheres are generated conditionally based on :yref:`withBoundaries` (not compatible with periodic boundaries)")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PartialSatClayEngine);

#endif //TwoPhaseFLOW
 
