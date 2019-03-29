/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              * 
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "FlowBoundingSphere.hpp"

#ifdef PARTIALSAT

namespace CGT {

//typedef CGT::_Tesselation<CGT::TriangulationTypes<PartialSatVertexInfo,PartialSatCellInfo> > PartialSatTesselation;

template<class _Tesselation>
class PartialSatLinSolv : public FlowBoundingSphereLinSolv<_Tesselation,FlowBoundingSphere<_Tesselation> >
{
public:
	typedef _Tesselation			Tesselation;
	typedef Network<Tesselation>		_N;
	DECLARE_TESSELATION_TYPES(Network<Tesselation>)
	typedef FlowBoundingSphereLinSolv<_Tesselation,FlowBoundingSphere<_Tesselation> > BaseFlowSolver;
	typedef typename BaseFlowSolver::ETriplet ETriplet; 
	
	///painfull, but we need that for templates inheritance...
	using _N::T; using _N::xMin; using _N::xMax; using _N::yMin; using _N::yMax; using _N::zMin; using _N::zMax; using _N::Rmoy; using _N::sectionArea; using _N::Height; using _N::vTotal; using _N::currentTes; using _N::debugOut; using _N::nOfSpheres; using _N::xMinId; using _N::xMaxId; using _N::yMinId; using _N::yMaxId; using _N::zMinId; using _N::zMaxId; using _N::boundsIds; using _N::cornerMin; using _N::cornerMax;  using _N::VSolidTot; using _N::Vtotalissimo; using _N::vPoral; using _N::sSolidTot; using _N::vPoralPorosity; using _N::vTotalPorosity; using _N::boundaries; using _N::idOffset; using _N::vtkInfiniteVertices; using _N::vtkInfiniteCells; using _N::num_particles; using _N::boundingCells; using _N::facetVertices; using _N::facetNFictious;
	//same for functions
	using _N::defineFictiousCells; using _N::addBoundingPlanes; using _N::boundary; ; using _N::tesselation; using _N::surfaceSolidThroatInPore;
	
	using BaseFlowSolver::noCache; using BaseFlowSolver::rAverage; using BaseFlowSolver::distanceCorrection; using BaseFlowSolver::minPermLength; using BaseFlowSolver::checkSphereFacetOverlap; using BaseFlowSolver::viscosity; using BaseFlowSolver::kFactor; using BaseFlowSolver::permeabilityMap; using BaseFlowSolver::maxKdivKmean; using BaseFlowSolver::clampKValues; using BaseFlowSolver::KOptFactor; using BaseFlowSolver::meanKStat; using BaseFlowSolver::fluidBulkModulus; using BaseFlowSolver::relax; using BaseFlowSolver::tolerance; using BaseFlowSolver::minKdivKmean; using BaseFlowSolver::resetRHS; using BaseFlowSolver::factorizeOnly; using BaseFlowSolver::perVertexUnitForce; using BaseFlowSolver::perVertexPressure; using BaseFlowSolver::computeAllCells; using BaseFlowSolver::setBlocked; using BaseFlowSolver::computeHydraulicRadius; using BaseFlowSolver::computeEquivalentRadius; using BaseFlowSolver::computeEffectiveRadius; 
	/// More members from LinSolv variant
	using BaseFlowSolver::areCellsOrdered; using BaseFlowSolver::T_nnz; using BaseFlowSolver::ncols; using BaseFlowSolver::T_cells; using BaseFlowSolver::T_index; using BaseFlowSolver::orderedCells; using BaseFlowSolver::isLinearSystemSet; using BaseFlowSolver::T_x; using BaseFlowSolver::T_b; using BaseFlowSolver::T_bv; using BaseFlowSolver::bodv; using BaseFlowSolver::xodv; using BaseFlowSolver::errorCode; using BaseFlowSolver::useSolver; using BaseFlowSolver::tripletList; using BaseFlowSolver::A; using BaseFlowSolver::gsP; using BaseFlowSolver::gsB; using BaseFlowSolver::fullAcolumns; using BaseFlowSolver::fullAvalues; using BaseFlowSolver::isFullLinearSystemGSSet; using BaseFlowSolver::gsdV; using BaseFlowSolver::multithread; using BaseFlowSolver::getCHOLMODPerfTimings; using BaseFlowSolver::thermalEngine; using BaseFlowSolver::fluidCp; using BaseFlowSolver::fluidRho; using BaseFlowSolver::reuseOrdering; using BaseFlowSolver::partialSatEngine; using BaseFlowSolver::factorExists; using BaseFlowSolver::start; using BaseFlowSolver::Achol; using BaseFlowSolver::com; using BaseFlowSolver::cholT; using BaseFlowSolver::L; using BaseFlowSolver::end; using BaseFlowSolver::updatedRHS; using BaseFlowSolver::ompThreads; using BaseFlowSolver::add_T_entry;
	
	vector<int> indices;//redirection vector containing the rank of cell so that T_cells[indices[cell->info().index]]=cell

	virtual ~PartialSatLinSolv();
	PartialSatLinSolv();

	///Linear system solve
	virtual int setLinearSystem(Real dt=0);
	virtual int setLinearSystemFullGS(Real dt=0);
	virtual void copyCellsToLin(Real dt=0);
	virtual void interpolate(Tesselation& Tes, Tesselation& NewTes);
	virtual void computeFacetForcesWithCache(bool onlyCache=false);
	virtual	void computePermeability();
	virtual	void gaussSeidel(Real dt=0);
};

} //namespace CGTF
#include "PartialSatLinSolv.ipp"
#endif //FLOW_ENGINE
