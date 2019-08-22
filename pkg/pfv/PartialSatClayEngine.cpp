/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              * 
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
// Experimental engine under development 
#ifdef PARTIALSAT
#include "PartialSatClayEngine.hpp"
#include <boost/range/algorithm_ext/erase.hpp>
#include<pkg/dem/CohesiveFrictionalContactLaw.hpp>
#include<pkg/common/Sphere.hpp>
#include<core/Body.hpp>

CREATE_LOGGER(PartialSatClayEngine);
YADE_PLUGIN((PartialSatClayEngine));

PartialSatClayEngine::~PartialSatClayEngine(){}

void PartialSatClayEngine::action(){
       if ( !isActivated ) return;
        timingDeltas->start();
        if (desiredPorosity != 0){
	    Real actualPorosity = Shop::getPorosityAlt();
	    volumeCorrection = desiredPorosity/actualPorosity;
        }	
	setPositionsBuffer(true);
	timingDeltas->checkpoint ( "Position buffer" );
        if (first) {
	  if (multithread) setPositionsBuffer(false);
	  buildTriangulation(pZero,*solver);
	  initializeVolumes(*solver);
	  backgroundSolver=solver;
	  backgroundCompleted=true;
          if (partialSatEngine) {
		initializeSaturations(*solver);
		if (particleSwelling && volumes) setOriginalParticleValues();
		if (debug) cout << "Particle swelling model active, original particle volumes set" << endl;
	  }
	}
	#ifdef YADE_OPENMP
	solver->ompThreads = ompThreads>0? ompThreads : omp_get_max_threads();
	#endif
        timingDeltas->checkpoint ( "Triangulating" );
	updateVolumes ( *solver );
        timingDeltas->checkpoint ( "Update_Volumes" );
	
        epsVolCumulative += epsVolMax;
	retriangulationLastIter++;
	if (!updateTriangulation) updateTriangulation = // If not already set true by another function of by the user, check conditions
		(defTolerance>0 && epsVolCumulative > defTolerance) || (meshUpdateInterval>0 && retriangulationLastIter>meshUpdateInterval);

	// remesh everytime a bond break occurs (for DFNFlow-JCFPM coupling)
	if (breakControlledRemesh) remeshForFreshlyBrokenBonds();

        ///compute flow and and forces here
	if (pressureForce){
	#ifdef LINSOLV
		permUpdateIters++;
		if (fixTriUpdatePermInt>0 && permUpdateIters >= fixTriUpdatePermInt) updateLinearSystem(*solver);
	#endif
		if (partialSatEngine) setCellsDSDP(*solver);
		solver->gaussSeidel(scene->dt);
		timingDeltas->checkpoint ( "Factorize + Solve" );
		if (!decoupleForces) solver->computeFacetForcesWithCache();
		if (partialSatEngine) {
			initializeSaturations(*solver);
			//updateSaturation(*solver);
			if (debug) cout << "finished initializing saturations" << endl;
		}
	}
	if (particleSwelling) {
		if (suction) collectParticleSuction(*solver);
		if (swelling) swellParticles();
	}
	if (debug) cout << "finished collecting suction and swelling " << endl;
	if (freeSwelling && crackModelActive) determineFracturePaths();
        timingDeltas->checkpoint ( "compute_Forces" );
        ///Application of vicscous forces
	scene->forces.sync();
	timingDeltas->checkpoint ( "forces.sync()" );
	if (!decoupleForces) computeViscousForces ( *solver );
	timingDeltas->checkpoint ( "viscous forces" );
	if (!decoupleForces) applyForces(*solver);
        timingDeltas->checkpoint ( "Applying Forces" );
	if (debug) cout << "finished cmputing forces and applying" <<endl;
	///End compute flow and forces
	#ifdef LINSOLV
	int sleeping = 0;
	if (multithread && !first) {
		while (updateTriangulation && !backgroundCompleted) { /*cout<<"sleeping..."<<sleeping++<<endl;*/
		  sleeping++;
		boost::this_thread::sleep(boost::posix_time::microseconds(1000));}
		if (debug && sleeping) cerr<<"sleeping..."<<sleeping<<endl;
		if (updateTriangulation || ((meshUpdateInterval>0  && ellapsedIter>(0.5*meshUpdateInterval)) && backgroundCompleted)) {
			if (debug) cerr<<"switch flow solver"<<endl;
			if (useSolver==0) LOG_ERROR("background calculations not available for Gauss-Seidel");
			if(!fluxChanged){
				if (fluidBulkModulus>0 || doInterpolate) solver->interpolate (solver->T[solver->currentTes], backgroundSolver->T[backgroundSolver->currentTes]);
				
				

			//Copy imposed pressures/flow from the old solver
				backgroundSolver->imposedP = vector<pair<CGT::Point,Real> >(solver->imposedP);
				backgroundSolver->imposedF = vector<pair<CGT::Point,Real> >(solver->imposedF);
				solver=backgroundSolver;
			} else {
				fluxChanged=false;
			}
			
			backgroundSolver = shared_ptr<FlowSolver> (new FlowSolver);
			if (metisForced) {backgroundSolver->eSolver.cholmod().nmethods=1; backgroundSolver->eSolver.cholmod().method[0].ordering=CHOLMOD_METIS;}
			backgroundSolver->imposedP = vector<pair<CGT::Point,Real> >(solver->imposedP);
			backgroundSolver->imposedF = vector<pair<CGT::Point,Real> >(solver->imposedF);
			if (debug) cerr<<"switched"<<endl;
			setPositionsBuffer(false);//set "parallel" buffer for background calculation 
			backgroundCompleted=false;
			retriangulationLastIter=ellapsedIter;
			updateTriangulation=false;
			epsVolCumulative=0;
			ellapsedIter=0;
			boost::thread workerThread(&PartialSatClayEngine::backgroundAction,this);
			workerThread.detach();
			if (debug) cerr<<"backgrounded"<<endl;
			initializeVolumes(*solver);
			computeViscousForces(*solver);
			if (debug) cerr<<"volumes initialized"<<endl;
		}
		else {
			if (debug && !backgroundCompleted) cerr<<"still computing solver in the background, ellapsedIter="<<ellapsedIter<<endl;
			ellapsedIter++;
		}
	} else
	#endif
	 {
	        if (updateTriangulation && !first) {
			buildTriangulation (pZero, *solver);
			initializeVolumes(*solver);
			computeViscousForces(*solver);
               		updateTriangulation = false;
			epsVolCumulative=0;
			retriangulationLastIter=0;
			ReTrg++;}
        }
        first=false;
        timingDeltas->checkpoint ( "triangulate + init volumes" );
}



double PartialSatClayEngine::diagonalSaturationContribution(CellHandle cell){
	double contr;
	if (cell->info().saturation<1.0) contr = cell->info().dsdp/(cell->info().invVoidVolume()*scene->dt);
	else contr = 0;
	return contr;
}

double PartialSatClayEngine::RHSSaturationContribution(CellHandle cell){
	double contr;
	if (cell->info().saturation<1.0) contr = cell->info().p() * cell->info().dsdp / (scene->dt * cell->info().invVoidVolume() );
	else contr = 0;
	return contr;
}

/////// Partial Sat Tools /////////

void PartialSatClayEngine::setCellsDSDP(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	#pragma omp parallel for
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
		const double deriv = dsdp(cell);
        	if (!isnan(deriv)) cell->info().dsdp = deriv;
		else cell->info().dsdp = 0;
	}
}	

double PartialSatClayEngine::dsdp(CellHandle& cell)
{
	// analytical derivative of van genuchten
	double pc = pAir - cell->info().p();
	double term1 = pow((pow((pc/Po),(1./1.-lmbda))+1.),(-lmbda-1.));
	double term2 = lmbda*pow((pc/Po),(1./(1.-lmbda)-1.));
	double term3 = Po*(1.-lmbda);
	return term1 * term2 / term3; 

}

void PartialSatClayEngine::initializeSaturations(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	#pragma omp parallel for
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
        	setSaturationFromPcS(cell);
	}
}

void PartialSatClayEngine::setSaturationFromPcS(CellHandle& cell)
{
	// using van genuchten	
	double pc = pAir - cell->info().p();
	double saturation;
	if (pc>=0) saturation = pow((1. + pow(pc/Po,1./(1.-lmbda))),-lmbda);
	else saturation = 1.;
	cell->info().saturation = saturation;
}

void PartialSatClayEngine::updateSaturation(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().Pcondition) continue;
		double sum=0;
		for (int j=0; j<4; j++){
			CellHandle neighborCell = cell->neighbor(j);
			sum += cell->info().kNorm()[j] * (cell->info().p() - neighborCell->info().p());
		}
        	cell->info().saturation = cell->info().saturation + scene->dt * cell->info().invVoidVolume() * sum;
	}
}

void PartialSatClayEngine::collectParticleSuction(FlowSolver& flow){
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	Tesselation& Tes = flow.T[flow.currentTes];
	const long size = Tes.cellHandles.size();
	#pragma omp parallel for
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().isGhost || cell->info().Pcondition || cell->info().isFictious) continue; // Do we need special cases for fictious cells?
		for (int v=0;v<4;v++){
			//if (cell->vertex(v)->info().isFictious) continue;
			const long int id = cell->vertex(v)->info().id();
			const shared_ptr<Body>& b =(*bodies)[id];
			if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
			auto* state = b->state.get();
			//if (cell->info().isExposed) state->suctionSum+= pAir; // use different pressure for exposed cracks?
			state->suctionSum += pAir - cell->info().p();
			state->incidentCells += 1;
		}
	}

}

void PartialSatClayEngine::setOriginalParticleValues(){
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
	#pragma omp parallel for
	for (long i=0; i<size; i++){
		const shared_ptr<Body>& b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		const double volume = 4./3. * M_PI*pow(sphere->radius,3.);
		auto* state = b->state.get();
		state->volumeOriginal = volume;
		state->radiiOriginal = sphere->radius;
	}
}


void PartialSatClayEngine::swellParticles(){
	const shared_ptr<BodyContainer>& bodies=scene->bodies;
	const long size=bodies->size();
	const double suction0 = pAir - pZero;
	totalVolChange = 0;
	#pragma omp parallel for
	for (long i=0; i<size; i++){
		const shared_ptr<Body>& b=(*bodies)[i];
		if (b->shape->getClassIndex()!=Sphere::getClassIndexStatic() || !b) continue;
		Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		//const double volume = 4./3. * M_PI*pow(sphere->radius,3);
		auto* state = b->state.get();
		state->suction = state->suctionSum/state->incidentCells;
		state->incidentCells = 0; // reset to 0 for next time step
		state->suctionSum = 0; //
		const double volStrain = betam/alpham * (exp(-alpham*state->suction) - exp(-alpham*suction0));
//		const double rOrig = pow(state->volumeOriginal * 3. / (4.*M_PI),1./3.);
	//
		const double vNew = state->volumeOriginal*(volStrain + 1.);
		const double rNew = pow(3.*vNew/(4.*M_PI),1./3.);
		totalVolChange += (pow(rNew,3.) - pow(sphere->radius,3.)) * 4./3.*M_PI;
		state->radiiChange = rNew - state->radiiOriginal;
		sphere->radius = rNew;
//		cout << "volStrain "<<volStrain<<" avgSuction "<<avgSuction<<" suction0 " <<suction0<<" rDel "<<rDel<<" rNew "<< rNew << " rOrig "<< rOrig << endl;
	}

}

//////// Post processing tools //////

void PartialSatClayEngine::saveUnsatVtk(const char* folder, bool withBoundaries)
{
	vector<int> allIds;//an ordered list of cell ids (from begin() to end(), for vtk table lookup), some ids will appear multiple times since boundary cells are splitted into multiple tetrahedra 
	vector<int> fictiousN;
	bool initNoCache=solver->noCache;
	solver->noCache=false;
	
	static unsigned int number = 0;
        char filename[250];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(filename,"%s/out_%d.vtk",folder,number++);
	basicVTKwritter vtkfile(0,0);
	solver->saveMesh(vtkfile,withBoundaries,allIds,fictiousN,filename);
	solver->noCache=initNoCache;
		
	vtkfile.begin_data("Pressure",CELL_DATA,SCALARS,FLOAT);
	for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().p());
	vtkfile.end_data();
		
	vtkfile.begin_data("fictious",CELL_DATA,SCALARS,INT);
	for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(fictiousN[kk]);
	vtkfile.end_data();

	vtkfile.begin_data("id",CELL_DATA,SCALARS,INT);
	for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(allIds[kk]);
	vtkfile.end_data();

	vtkfile.begin_data("fracturedCells",CELL_DATA,SCALARS,INT);
	for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().crack);
	vtkfile.end_data();


	vtkfile.begin_data("Permeability",CELL_DATA,SCALARS,FLOAT);
	for (unsigned kk=0; kk<allIds.size(); kk++) {
		std::vector<double> perm = solver->tesselation().cellHandles[allIds[kk]]->info().kNorm();
		double permSum = 0;
		for (unsigned int i=0;i<perm.size();i++) permSum+=perm[i];
		vtkfile.write_data(permSum);
	}
	vtkfile.end_data();
	
	#define SAVE_CELL_INFO(INFO) vtkfile.begin_data(#INFO,CELL_DATA,SCALARS,FLOAT); for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().INFO); vtkfile.end_data();
	SAVE_CELL_INFO(saturation)
	SAVE_CELL_INFO(Pcondition)
	SAVE_CELL_INFO(isExposed)
	//SAVE_CELL_INFO(porosity)
}

void PartialSatClayEngine::initSolver ( FlowSolver& flow )
{
       	flow.Vtotalissimo=0; flow.VSolidTot=0; flow.vPoral=0; flow.sSolidTot=0;
        flow.slipBoundary=slipBoundary;
        flow.kFactor = permeabilityFactor;
        flow.debugOut = debug;
        flow.useSolver = useSolver;
	#ifdef CHOLMOD_LIBS
	flow.numSolveThreads = numSolveThreads;
	flow.numFactorizeThreads = numFactorizeThreads;
	#endif
	flow.factorizeOnly = false;
	flow.meanKStat = meanKStat;
        flow.viscosity = viscosity;
        flow.tolerance=tolerance;
        flow.relax=relax;
        flow.clampKValues = clampKValues;
	flow.maxKdivKmean = maxKdivKmean;
	flow.minKdivKmean = minKdivKmean;
        flow.meanKStat = meanKStat;
        flow.permeabilityMap = permeabilityMap;
        flow.fluidBulkModulus = fluidBulkModulus;
	flow.fluidRho = fluidRho;
	flow.fluidCp = fluidCp;
	flow.thermalEngine = thermalEngine;
	flow.multithread = multithread;
	flow.getCHOLMODPerfTimings=getCHOLMODPerfTimings;
//         flow.tesselation().Clear();
        flow.tesselation().maxId=-1;
	flow.blockedCells.clear();
	flow.sphericalVertexAreaCalculated=false;
        flow.xMin = 1000.0, flow.xMax = -10000.0, flow.yMin = 1000.0, flow.yMax = -10000.0, flow.zMin = 1000.0, flow.zMax = -10000.0;
	flow.partialSatEngine = partialSatEngine;
	flow.pAir = pAir;
	flow.freeSwelling = freeSwelling;
	flow.matricSuctionRatio = matricSuctionRatio;
	flow.nUnsatPerm = nUnsatPerm;
	flow.SrM = SrM;
	flow.SsM = SsM;
}

void PartialSatClayEngine::buildTriangulation ( double pZero, Solver& flow )
{
 	if (first) flow.currentTes=0;
        else {  flow.currentTes=!flow.currentTes; if (debug) cout << "--------RETRIANGULATION-----------" << endl;}
	flow.resetNetwork();
	initSolver(flow);

        addBoundary ( flow );
        triangulate ( flow );
        if ( debug ) cout << endl << "Tesselating------" << endl << endl;
        flow.tesselation().compute();
        flow.defineFictiousCells();
	// For faster loops on cells define this vector
	flow.tesselation().cellHandles.clear();
	flow.tesselation().cellHandles.reserve(flow.tesselation().Triangulation().number_of_finite_cells());
	FiniteCellsIterator cell_end = flow.tesselation().Triangulation().finite_cells_end();
	int k=0;
	for ( FiniteCellsIterator cell = flow.tesselation().Triangulation().finite_cells_begin(); cell != cell_end; cell++ ){
		flow.tesselation().cellHandles.push_back(cell);
		cell->info().id=k++;}//define unique numbering now, corresponds to position in cellHandles
        flow.displayStatistics ();
	if(!blockHook.empty()){ LOG_INFO("Running blockHook: "<<blockHook); pyRunString(blockHook); }
        flow.computePermeability();

	if (multithread && fluidBulkModulus>0) initializeVolumes(flow);  // needed for multithreaded compressible flow (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1687355)
//	if (crackModelActive) trickPermeability(&flow); //This virtual function does nothing yet, derived class may overload it to make permeability different (see DFN engine)
        porosity = flow.vPoralPorosity/flow.vTotalPorosity;

        boundaryConditions ( flow );
        flow.initializePressure ( pZero );
	if (thermalEngine) {
		//initializeVolumes(flow);
		thermalBoundaryConditions ( flow ); 
		flow.initializeTemperatures ( tZero );
		flow.sphericalVertexAreaCalculated=false;
	}

	
        if ( !first && !multithread && (useSolver==0 || fluidBulkModulus>0 || doInterpolate || thermalEngine || partialSatEngine)) flow.interpolate ( flow.T[!flow.currentTes], flow.tesselation() );
	if (crackModelActive) trickPermeability(&flow);
        if ( waveAction ) flow.applySinusoidalPressure ( flow.tesselation().Triangulation(), sineMagnitude, sineAverage, 30 );
	else if (boundaryPressure.size()!=0) flow.applyUserDefinedPressure ( flow.tesselation().Triangulation(), boundaryXPos , boundaryPressure);
        if (normalLubrication || shearLubrication || viscousShear) flow.computeEdgesSurfaces();
}


void PartialSatClayEngine::initializeVolumes (FlowSolver& flow )
{
	typedef typename Solver::FiniteVerticesIterator FiniteVerticesIterator;
	
	FiniteVerticesIterator vertices_end = flow.tesselation().Triangulation().finite_vertices_end();
	CGT::CVector Zero(0,0,0);
	for (FiniteVerticesIterator V_it = flow.tesselation().Triangulation().finite_vertices_begin(); V_it!= vertices_end; V_it++) V_it->info().forces=Zero;

	FOREACH(CellHandle& cell, flow.tesselation().cellHandles)
	{
		switch ( cell->info().fictious() )
		{
			case ( 0 ) : cell->info().volume() = volumeCell ( cell ); break;
			case ( 1 ) : cell->info().volume() = volumeCellSingleFictious ( cell ); break;
			case ( 2 ) : cell->info().volume() = volumeCellDoubleFictious ( cell ); break;
			case ( 3 ) : cell->info().volume() = volumeCellTripleFictious ( cell ); break;
			default: break; 
		}
	if (flow.fluidBulkModulus>0 || thermalEngine || iniVoidVolumes) {
		cell->info().invVoidVolume() = 1 / ( std::abs(cell->info().volume()) - volumeCorrection*flow.volumeSolidPore(cell) );
	} else if (partialSatEngine) { 
		cell->info().invVoidVolume() = 1 / std::abs(cell->info().volume()); 
	}
	}
	if (debug) cout << "Volumes initialised." << endl;
}



/////// Discrete Fracture Network Functionality ////////


void PartialSatClayEngine::interpolateCrack(Tesselation& Tes,Tesselation& NewTes){
        RTriangulation& Tri = Tes.Triangulation();
	//RTriangulation& newTri = NewTes.Triangulation();
	//FiniteCellsIterator cellEnd = newTri.finite_cells_end();
	#ifdef YADE_OPENMP
    	const long size = NewTes.cellHandles.size();
	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& newCell = NewTes.cellHandles[i];
        #else
	FOREACH(CellHandle& newCell, NewTes.cellHandles){
        #endif
		if (newCell->info().isGhost) continue;
		CVector center ( 0,0,0 );
		if (newCell->info().fictious()==0) for ( int k=0;k<4;k++ ) center= center + 0.25* (Tes.vertex(newCell->vertex(k)->info().id())->point().point()-CGAL::ORIGIN);
		else {
			Real boundPos=0; int coord=0;
			for ( int k=0;k<4;k++ ){
				if (!newCell->vertex (k)->info().isFictious) center= center+(1./(4.-newCell->info().fictious()))*(Tes.vertex(newCell->vertex(k)->info().id())->point().point()-CGAL::ORIGIN);
			}
			for ( int k=0;k<4;k++ ) {
				if (newCell->vertex (k)->info().isFictious) {
					coord=solver->boundary (newCell->vertex(k)->info().id()).coordinate;
					boundPos=solver->boundary (newCell->vertex(k)->info().id()).p[coord];
					center=CVector(coord==0?boundPos:center[0],coord==1?boundPos:center[1],coord==2?boundPos:center[2]);
				}
			}
		}
		CellHandle oldCell = Tri.locate(CGT::Sphere(center[0],center[1],center[2]));
		newCell->info().crack = oldCell->info().crack;
//		For later commit newCell->info().fractureTip = oldCell->info().fractureTip;
//		For later commit newCell->info().cellHalfWidth = oldCell->info().cellHalfWidth;

		/// compute leakoff rate by summing the flow through facets abutting non-cracked neighbors
//		if (oldCell->info().crack && !oldCell->info().fictious()){
//			Real facetFlowRate=0;
//			facetFlowRate -= oldCell->info().dv();
//			for (int k=0; k<4;k++) {
//				if (!oldCell->neighbor(k)->info().crack){
//					facetFlowRate = oldCell->info().kNorm()[k]*(oldCell->info().shiftedP()-oldCell->neighbor(k)->info().shiftedP());
//					leakOffRate += facetFlowRate;
//				}
//			}
//		}
	}
    }

void PartialSatClayEngine::computeFracturePerm(RTriangulation::Facet_circulator& facet,Real aperture,RTriangulation::Finite_edges_iterator& ed_it)
{
	const RTriangulation::Facet& currentFacet = *facet; /// seems verbose but facet->first was declaring a junk cell and crashing program (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1666339)
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	const CellHandle& cell1 = currentFacet.first;
	const CellHandle& cell2 = currentFacet.first->neighbor(facet->second);
	if (Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr<<"Infinite cell found in trickPermeability, should be handled somehow, maybe"<<endl;
	
// 	/// ROBERT
// 	if (cell1->info().count()[currentFacet.second] < 3){
// 		cell1->info().count()[currentFacet.second] += 1;
// 		//cell1->info().aperture()[currentFacet.second] += aperture // used with avgAperture below if desired
// 	}
// 	if (!edgeOnJoint && cell1->info().count()[currentFacet.second] < facetEdgeBreakThreshold) return; // only allow facets with 2 or 3 broken edges to be tricked
	// Need to decide aperture criteria! Otherwise it selects the most recent broken edge aperture for facet's perm calc. Here is one possibility: 
	//Real avgAperture = cell1->info().aperture()[currentFacet.second]/Real(cell1->info().count()[currentFacet.second]);
	double fracturePerm = apertureFactor*pow(aperture,3.)/(12.*viscosity);
	cell1->info().kNorm()[currentFacet.second] += fracturePerm; //
	cell2->info().kNorm()[Tri.mirror_index(cell1,currentFacet.second)] += fracturePerm;
	/// For vtk recorder:
	cell1->info().crack=1;
	cell2->info().crack=1;
	//cout << "crack set to true in pore"<<endl;
	cell2->info().blocked = cell1->info().blocked = cell2->info().Pcondition = cell1->info().Pcondition = false; /// those ones will be included in the flow problem
	Point& CellCentre1 = cell1->info(); /// Trying to get fracture's surface 
	Point& CellCentre2 = cell2->info(); /// Trying to get fracture's surface 
	CVector networkFractureLength = CellCentre1 - CellCentre2; /// Trying to get fracture's surface 
	double networkFractureDistance = sqrt(networkFractureLength.squared_length()); /// Trying to get fracture's surface 
	Real networkFractureArea = pow(networkFractureDistance,2);  /// Trying to get fracture's surface 
	totalFractureArea += networkFractureArea; /// Trying to get fracture's surface 
// 	cout <<" ------------------ The total surface area up to here is --------------------" << totalFractureArea << endl;
// 	printFractureTotalArea = totalFractureArea; /// Trying to get fracture's surface 
	if (calcCrackArea) {
			CVector edge = ed_it->first->vertex(ed_it->second)->point().point() - ed_it->first->vertex(ed_it->third)->point().point();
			CVector unitV = edge*(1./sqrt(edge.squared_length()));
			Point p3 = ed_it->first->vertex(ed_it->third)->point().point() + unitV*(cell1->info() - ed_it->first->vertex(ed_it->third)->point().point())*unitV;
			Real halfCrackArea = 0.25*sqrt(std::abs(cross_product(CellCentre1-p3,CellCentre2-p3).squared_length()));//
			cell1->info().crackArea += halfCrackArea;
			cell2->info().crackArea += halfCrackArea;
			crackArea += 2*halfCrackArea;
		}
}

void PartialSatClayEngine::circulateFacets(RTriangulation::Finite_edges_iterator& edge,Real aperture)
{
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0 = facet1++;
	computeFracturePerm(facet0,aperture,edge);
	while (facet1!=facet0) {computeFracturePerm(facet1,aperture,edge); facet1++;}
	/// Needs the fracture surface for this edge?
	// 	double edgeArea = solver->T[solver->currentTes].computeVFacetArea(edge); cout<<"edge area="<<edgeArea<<endl;
}

void PartialSatClayEngine::trickPermeability(Solver* flow)
{
	leakOffRate=0;	
	const RTriangulation& Tri = flow->T[solver->currentTes].Triangulation();
//	if (!first) interpolateCrack(solver->T[solver->currentTes],flow->T[flow->currentTes]);
	const CohFrictPhys* cohfrictphys;
	const shared_ptr<InteractionContainer> interactions = scene->interactions;
	int numberOfCrackedOrJointedInteractions=0; 
	Real SumOfApertures=0.; 
	averageAperture=0;
	maxAperture=0;
	crackArea=0;
	//Real totalFractureArea=0; /// Trying to get fracture's surface
// 	const shared_ptr<IGeom>& ig;
// 	const ScGeom* geom; // = static_cast<ScGeom*>(ig.get());
	FiniteEdgesIterator edge = Tri.finite_edges_begin();
	for( ; edge != Tri.finite_edges_end(); ++edge) {
	  
		const VertexInfo& vi1 = (edge->first)->vertex(edge->second)->info();
		const VertexInfo& vi2 = (edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction = interactions->find( vi1.id(),vi2.id() );
		
		if (interaction && interaction->isReal()) {
			if (edge->first->info().isFictious) continue; /// avoid trick permeability for fictitious
			cohfrictphys = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
			
			if (cohfrictphys->isBroken || allCellsFractured) {
				numberOfCrackedOrJointedInteractions+=1;
				//if (cohfrictphys->crackAperture<=0.) continue;
				Real aperture = (cohfrictphys->crackAperture <= residualAperture)? residualAperture : cohfrictphys->crackAperture;

				if (aperture > maxAperture) maxAperture = aperture; 
				SumOfApertures += aperture;
				
				circulateFacets(edge,aperture);
			};
		}
	}
	averageAperture = SumOfApertures/numberOfCrackedOrJointedInteractions; /// DEBUG
// 	cout << " Average aperture in joint ( -D ) = " << AverageAperture << endl; /// DEBUG
}

void PartialSatClayEngine::determineFracturePaths()
{
    RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
    FiniteCellsIterator cellEnd = tri.finite_cells_end();
    for ( FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++ ) {
      if(cell->info().Pcondition) continue;
        cell->info().isExposed = false;
    }
    
    for (int i=0;i<6;i++){
    	for (FlowSolver::VCellIterator it = solver->boundingCells[i].begin(); it != solver->boundingCells[i].end(); it++) {
      		if ((*it)==NULL) continue;
        	exposureRecursion(*it);
    	}
    }
}

void PartialSatClayEngine::exposureRecursion(CellHandle cell)
{
    for (int facet = 0; facet < 4; facet ++) {
        CellHandle nCell = cell->neighbor(facet);
        if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell)) continue;
        if (nCell->info().Pcondition) continue;
//         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
        if (!nCell->info().crack) continue;
       	if (nCell->info().isExposed) continue; // another recursion already found it
        nCell->info().isExposed=true;
	//nCell->info().p()=pAir;	
        exposureRecursion(nCell);
    }
}

bool PartialSatClayEngine::findInscribedRadiusAndLocation(CellHandle& cell, std::vector<double>& coordAndRad)
{
	//cout << "using least sq to find inscribed radius " << endl;
	const float prec = 1e-5;
	Eigen::MatrixXd A(4,3);
	Eigen::Vector4d b;
	Eigen::Vector3d x;
	Eigen::Vector4d r;
	//std::vector<double> r(4);
	//std::vector<double> coordAndRad(4); 
	CVector baryCenter ( 0,0,0 ); // use cell barycenter as initial guess
	for ( int k=0;k<4;k++ ){
		baryCenter= baryCenter + 0.25* (cell->vertex(k)->point().point()-CGAL::ORIGIN);
		if (cell->vertex(k)->info().isFictious) return 0;
	}
	double xo,yo,zo;
	int count=0;
	double rMean;
	xo=baryCenter[0];yo=baryCenter[1];zo=baryCenter[2];
	bool finished = false;
	while (finished==false) { 
		count += 1;
		if (count>1000) {
			cerr<<"too many iterations during sphere inscription, quitting"<<endl;
			return 0;
		}
		// build A matrix (and part of b)
		for (int i=0;i<4;i++){
			double xi,yi,zi;
			xi=cell->vertex(i)->point().x();yi=cell->vertex(i)->point().y();zi=cell->vertex(i)->point().z();
			A(i,0) = xo - cell->vertex(i)->point().x(); 
			A(i,1) = yo - cell->vertex(i)->point().y(); 
			A(i,2) = zo - cell->vertex(i)->point().z();
			const double sqrdD = pow(xo - xi,2) + pow(yo - yi,2) + pow(zo - zi,2);
			r(i) = sqrt(sqrdD) - sqrt(cell->vertex(i)->point().weight());
		}	
		rMean = r.sum()/4.;
	
		// build b 
		for (int i=0;i<4;i++){
			double xi,yi,zi;
			xi=cell->vertex(i)->point().x();yi=cell->vertex(i)->point().y();zi=cell->vertex(i)->point().z();
			b(i) =( pow(rMean+sqrt(cell->vertex(i)->point().weight()),2.) - (pow(xo-xi,2.) + pow(yo-yi,2.) + pow(zo-zi,2.)) )/2.;
		}
	
		// use least squares (normal equation) to minimize residuals
		x = (A.transpose() * A).ldlt().solve(A.transpose()*b);
		// if the values are greater than precision, update the guess and repeat
		if (abs(x(0))>prec || abs(x(1))>prec || abs(x(2))>prec){
			xo += x(0);
			yo += x(1);
			zo += x(2);
		} else {
			coordAndRad[0]=xo+x(0);coordAndRad[1]=yo+x(1);coordAndRad[2]=zo+x(2);coordAndRad[3]=rMean;
			if (rMean > sqrt(cell->vertex(0)->point().weight())) return 0; // inscribed sphere might be excessively large if it is in a flat boundary cell
			finished = true;
		}
	
	} // end while finished == false
	
	return 1;

}

void PartialSatClayEngine::insertMicroPores(const float fracMicroPore)
{
	cout << "Inserting micro pores in " << fracMicroPore << " perc. of existing pores " << endl;
	Eigen::MatrixXd M(6,6);
	//if (!solver->T[solver->currentTes]){cerr << "No triangulation, not inserting micropores" << endl; return}
	Tesselation& Tes = solver->T[solver->currentTes];
	//cout << "Tes set" << endl;
	const long size = Tes.cellHandles.size();
	std::vector<bool> visited(size);
	std::vector<int> poreIndices(int(ceil(fracMicroPore*size)));
	bool numFound;
	// randomly select the pore indices that we will turn into micro pores
	#pragma omp parallel for
    	for (unsigned int i=0; i<poreIndices.size(); i++){
		numFound=false;
		while (!numFound) { 
			const double num = rand() % size;  // + 1?
			if ( !visited[num] && !Tes.cellHandles[num]->info().isFictious){
				visited[num]=true;
				poreIndices[i]=num;
				numFound=true;
			}
		}
	}
	cout << "find inscribed sphere radius" << endl;

	// find inscribed sphere radius in selected pores and add body 
	// FIXME How do we deal with inscribed spheres that might be overlapping after inscription?
	std::vector<double> coordsAndRad(4);
	//#pragma omp parallel for
      	for (unsigned int i=0; i<poreIndices.size(); i++){
		const int idx = poreIndices[i];
		CellHandle& cell = Tes.cellHandles[idx];
		for (int j=0;j<4;j++) if(cell->neighbor(j)->info().isFictious) continue; // avoid inscribing spheres in very flat exterior cells (FIXME: we can avoid this by using a proper alpha shape)
		//if (cell->info().Pcondition) continue;
		bool inscribed = findInscribedRadiusAndLocation(cell,coordsAndRad);
		if (!inscribed) continue; // sphere couldn't be inscribed, continue loop
		bool contained = checkSphereContainedInTet(cell,coordsAndRad);
		if (!contained) continue;
		//cout << "converting to Vector3r" << endl;
		Vector3r coords; coords[0]=Real(coordsAndRad[0]); coords[1]=Real(coordsAndRad[1]); coords[2]=Real(coordsAndRad[2]);
		const double radius = coordsAndRad[3];
		//cout << "adding body" << endl;
		shared_ptr<Body> body;
		createSphere(body,coords,radius,false,true);
		scene->bodies->insert(body);
	}
}

//bool PartialSatClayEngine::checkSphereContainedInTet(CellHandle& cell,std::vector<double>& coordsAndRad)
//{
//	Eigen::Vector3d inscSphere(coordsAndRad[0],coordsAndRad[1],coordsAndRad[2]);
//	Eigen::Vector3d cellLoc(cell->info()[0],cell->info()[1],cell->info()[2]);
//	double radius = coordsAndRad[3];
//	 for ( int i=0; i<4; i++ ) {
//		Eigen::Vector3d neighborCellLoc(cell->neighbor(i)->info()[0],cell->neighbor(i)->info()[1],cell->neighbor(i)->info()[2]);
//		Eigen::Vector3d vertexLoc(cell->vertex(facetVertices[i][0])->point().x(),cell->vertex(facetVertices[i][0])->point().y(),cell->vertex(facetVertices[i][0])->point().z());

//		Eigen::Vector3d Surfk = cellLoc-neighborCellLoc;
//		Eigen::Vector3d SurfkNormed = Surfk.normalized();
//		Eigen::Vector3d branch = vertexLoc - inscSphere;
//		double distToFacet = branch.dot(SurfkNormed);
//		if (distToFacet<0){
//			cerr<< "sphere center outside tet, skipping insertion"<<endl;
//			return 0;
//		} else if (distToFacet<radius) {
//			cerr << "inscribed sphere too large for tetrahedral, decreasing size from "<< radius <<" to "<<distToFacet<<endl;
//			coordsAndRad[3] = distToFacet;
//			radius = distToFacet;
//		}
//	}
//	return 1;		
//}

bool PartialSatClayEngine::checkSphereContainedInTet(CellHandle& cell,std::vector<double>& coordsAndRad)
{
	Eigen::Vector3d inscSphere(coordsAndRad[0],coordsAndRad[1],coordsAndRad[2]);
	const double origRadius = coordsAndRad[3];
	//Eigen::Vector3d neighborCellLoc(cell->neighbor(i)->info()[0],cell->neighbor(i)->info()[1],cell->neighbor(i)->info()[2]);
//	Eigen::MatrixXd A(3,4);
//	Eigen::Vector4d x;
//	Eigen::Vector3d bvec(0,0,0);
//	double a,b,c,d;
	double radius = coordsAndRad[3];
	 for ( int i=0; i<4; i++ ) {

		// using same logic as above but more explicit 
		Eigen::Vector3d nhat(cell->info().facetSurfaces[i][0],cell->info().facetSurfaces[i][1],cell->info().facetSurfaces[i][2]);
		nhat = nhat/sqrt(cell->info().facetSurfaces[i].squared_length());
		Eigen::Vector3d xi(cell->vertex(facetVertices[i][0])->point().x(),cell->vertex(facetVertices[i][0])->point().y(),cell->vertex(facetVertices[i][0])->point().z());
		double distToFacet = nhat.dot(inscSphere-xi);
		double exampleScale = sqrt(cell->vertex(facetVertices[i][0])->point().weight());
		double scale = exampleScale*minMicroRadFrac;
		// even more explicit, creating plane out of 3 verticies!
//		for (int j=0;j<3;j++){
//			A(j,0) = cell->vertex(facetVertices[i][j])->point().x();
//			A(j,1) = cell->vertex(facetVertices[i][j])->point().y();
//			A(j,2) = cell->vertex(facetVertices[i][j])->point().z();
//			A(j,3) = 1;
//		}
		if(!(distToFacet>=scale)) {
			cout<<"minimum radius size doesn't fit,in tet skipping"<<endl;
			return 0;
		}
//		x = A.colPivHouseholderQr().solve(bvec);
//		a=x(0);b=x(1);c=x(2);d=x(3);
//		double sqrtSum = sqrt(a*a+b*b+c*c);
//		double distToFacet = (a*coordsAndRad[0]+b*coordsAndRad[1]+c*coordsAndRad[2]+d)/sqrtSum; 

		if (distToFacet<0){
			cerr<< "sphere center outside tet, skipping insertion"<<endl;
			return 0;
		} else if (distToFacet<radius) {
			cerr << "inscribed sphere too large for tetrahedral, decreasing size from "<< radius <<" to "<<distToFacet<<endl;
			coordsAndRad[3] = distToFacet;//*(1.-minMicroRadFrac);
			radius = distToFacet;//*(1.-minMicroRadFrac);
		} //else {
			//cerr << "inscribed sphere too small, skipping insertion, btw rad*minMicro= " << exampleScale*minMicroRadFrac << " while dist to facet = " << distToFacet << " and the logic " << (distToFacet>=scale) << endl;
		//	return 0;
		//}
	}
	return 1;
}

void PartialSatClayEngine::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool big, bool dynamic )
{
	body = shared_ptr<Body>(new Body); body->groupMask=2;
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<Sphere> iSphere(new Sphere);
	body->state->blockedDOFs=State::DOF_NONE;
	const double volume = 4./3. * M_PI*pow(radius,3.);
	body->state->mass	= volume*microStructureRho;
	//body->state->inertia	= Vector3r(2.0/5.0*body->state->mass*radius*radius,
	//			2.0/5.0*body->state->mass*radius*radius,
	//			2.0/5.0*body->state->mass*radius*radius);
	body->state->pos=position;
	body->state->volumeOriginal = volume;
	body->state->radiiOriginal = radius;
	shared_ptr<CohFrictMat> mat(new CohFrictMat);
	mat->young		= microStructureE;
	mat->poisson		= microStructureNu;
	mat->frictionAngle	= microStructurePhi * Mathr::PI/180.0; //compactionFrictionDeg * Mathr::PI/180.0;
	mat->normalCohesion = mat->shearCohesion = microStructureAdh;	
	aabb->color		= Vector3r(0,1,0);
	iSphere->radius		= radius;
	//iSphere->color	= Vector3r(0.4,0.1,0.1);
	//iSphere->color           = Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom());
	//iSphere->color.normalize();
	body->shape	= iSphere;
	body->bound	= aabb;
	body->material	= mat;
}

//double PartialSatClayEngine::dsdp(CellHandle cell, double pw)
//{
//  
//    
//      if(pw == 0){std::cout << endl << "Error! water pressure is zero, while computing capillary pressure ... cellId= "<< cell->info().id;}
//      double exp = std::exp(-1*getKappa(cell->info().numberFacets) * cell->info().saturation);
//      double dsdp2 = (1.0 / cell->info().thresholdPressure) * std::pow((1.0 - exp),2.0) / ( getKappa(cell->info().numberFacets) * exp);
////       if(std::abs(dsdp2) > 1e10){ std::cerr << "Huge dsdp! : "<< dsdp2 << " " << exp << " "<< cell->info().thresholdPressure << " " << getKappa(cell->info().numberFacets);}

//    
////     double dsdp2 = (3.0 * getConstantC3(cell) - 2.0 * getConstantC4(cell) * pw) / std::pow(pw,4);
//    if(dsdp2 != dsdp2){std::cerr<<endl << "Error! sat in dsdp is nan: " << cell->info().saturation << " kappa:" <<getKappa(cell->info().numberFacets) << " exp: " << exp <<   " mergedVolume=" << cell->info().mergedVolume  << " pthreshold=" << cell->info().thresholdPressure;}
//    if(dsdp2 < 0.0){std::cerr<<endl<< "Error! dsdp is negative!" << dsdp2; dsdp2 = 0.0;}
////     if(dsdp2 >1e6){std::cerr<<endl<< "Error! dsdp is huge!" << dsdp2; dsdp2 = 1e6;}
//    return dsdp2;
//}

//double PartialSatClayEngine::poreSaturationFromPcS(CellHandle cell,double pw)
//{
//    //Using equation: Pc = 2*surfaceTension / (Chi * PoreBodyVolume^(1/3) * (1-exp(-kappa * S)))
//  double s = truncationPrecision;
//     if(-1*pw > cell->info().thresholdPressure){
//       s = std::log(1.0 + cell->info().thresholdPressure / pw) / (-1.0 * getKappa(cell->info().numberFacets));
//     }
//     if(-1*pw == cell->info().thresholdPressure){
//	s = cell->info().thresholdSaturation;
//     }
//     if(-1*pw < cell->info().thresholdPressure){
//	if(!remesh && !firstDynTPF){std::cerr<<endl<< "Error! Requesting saturation while capillary pressure is below threshold value? " << pw << " " << cell->info().thresholdPressure;}
//	s = cell->info().thresholdSaturation;
//     }

//     if(s > 1.0 || s < 0.0){std::cout << "Error, saturation from Pc(S) curve is not correct: "<< s << " "<<  cell->info().poreId << " log:" << std::log(1.0 + cell->info().thresholdPressure / pw)<< " " << (-1.0 * getKappa(cell->info().numberFacets)) << " pw=" << pw << " " << cell->info().thresholdPressure; s = 1.0;}
//     if(s != s){std::cerr<<endl << "Error! sat in PcS is nan: " << s << "  " << pw << " " << getConstantC4(cell) << " " << getConstantC3(cell) << " mergedVolume=" << cell->info().mergedVolume << " pthreshold=" << cell->info().thresholdPressure;}
//  return s;
//}


//double PartialSatClayEngine::porePressureFromPcS(CellHandle cell,double saturation)
//{
//   
//  double pw = -1.0 * cell->info().thresholdPressure / (1.0 - std::exp(-1*getKappa(cell->info().numberFacets) * cell->info().saturation));
//  if(std::exp(-1*getKappa(cell->info().numberFacets) * cell->info().saturation) == 1.0){std::cerr << endl << "Error! pw = -inf!"  << cell->info().saturation;}
//  if(pw > 0){
//    std::cout << "Pw is above 0! - error: "<< pw << " id=" << cell->info().poreId << " pthr=" << cell->info().thresholdPressure << " sat:" << cell->info().saturation << " kappa: " << getKappa(cell->info().numberFacets) << " " << (1.0 - std::exp(-1*getKappa(cell->info().numberFacets) * cell->info().saturation));
//    pw = -1 * cell->info().thresholdPressure;
//  }
//  if(pw != pw){std::cout << "Non existing capillary pressure!";}
////   if(pw < 100 * waterBoundaryPressure){std::cout << "huge PC!" << pw << " saturation=" << cell->info().saturation << " hasIFace=" << cell->info().hasInterface << " NWRES=" << cell->info().isNWRes; pw = 100 * waterBoundaryPressure;}
//  return pw; 
//}


#endif //PartialSat


