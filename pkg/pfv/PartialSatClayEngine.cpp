/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              * 
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef PARTIALSAT
#include "PartialSatClayEngine.hpp"
#include <boost/range/algorithm_ext/erase.hpp>

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
          if (partialSatEngine) initializeSaturations(*solver);
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
		solver->computeFacetForcesWithCache();
		if (partialSatEngine) {
			initializeSaturations(*solver);
			//updateSaturation(*solver);
		}
	}
        timingDeltas->checkpoint ( "compute_Forces" );
        ///Application of vicscous forces
        if (!decoupleForces) scene->forces.sync();
	timingDeltas->checkpoint ( "forces.sync()" );
	if (!decoupleForces) computeViscousForces ( *solver );
	timingDeltas->checkpoint ( "viscous forces" );
	if (!decoupleForces) applyForces(*solver);
        timingDeltas->checkpoint ( "Applying Forces" );
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
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& cell = Tes.cellHandles[i];
        	cell->info().dsdp = dsdp(cell);
	}
}	

double PartialSatClayEngine::dsdp(CellHandle& cell)
{
	// analytical derivative of van genuchten
	double pc = pAir - cell->info().p();
	double term1 = -lmbda*pow((pow((pc/Po),(1./1.-lmbda)) + 1.),-lmbda-1.);
	double term2 = pow(pc/Po,1./(1.-lmbda)-1.);
	double term3 = Po*(1.-lmbda);
	return term1 * term2 / term3; 

}

void PartialSatClayEngine::initializeSaturations(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
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
	
	#define SAVE_CELL_INFO(INFO) vtkfile.begin_data(#INFO,CELL_DATA,SCALARS,FLOAT); for (unsigned kk=0; kk<allIds.size(); kk++) vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().INFO); vtkfile.end_data();
	SAVE_CELL_INFO(saturation)
	SAVE_CELL_INFO(Pcondition)
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
}

//void PartialSatClayEngine::buildTriangulation ( double pZero, Solver& flow )
//{
// 	if (first) flow.currentTes=0;
//        else {  flow.currentTes=!flow.currentTes; if (debug) cout << "--------RETRIANGULATION-----------" << endl;}
//	flow.resetNetwork();
//	initSolver(flow);

//        addBoundary ( flow );
//        triangulate ( flow );
//        if ( debug ) cout << endl << "Tesselating------" << endl << endl;
//        flow.tesselation().compute();
//        flow.defineFictiousCells();
//	// For faster loops on cells define this vector
//	flow.tesselation().cellHandles.clear();
//	flow.tesselation().cellHandles.reserve(flow.tesselation().Triangulation().number_of_finite_cells());
//	FiniteCellsIterator cell_end = flow.tesselation().Triangulation().finite_cells_end();
//	int k=0;
//	for ( FiniteCellsIterator cell = flow.tesselation().Triangulation().finite_cells_begin(); cell != cell_end; cell++ ){
//		flow.tesselation().cellHandles.push_back(cell);
//		cell->info().id=k++;}//define unique numbering now, corresponds to position in cellHandles
//        flow.displayStatistics ();
//	if(!blockHook.empty()){ LOG_INFO("Running blockHook: "<<blockHook); pyRunString(blockHook); }
//        flow.computePermeability();

//	if (multithread && fluidBulkModulus>0) initializeVolumes(flow);  // needed for multithreaded compressible flow (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1687355)
//	trickPermeability(&flow); //This virtual function does nothing yet, derived class may overload it to make permeability different (see DFN engine)
//        porosity = flow.vPoralPorosity/flow.vTotalPorosity;

//        boundaryConditions ( flow );
//        flow.initializePressure ( pZero );
//	if (thermalEngine) {
//		//initializeVolumes(flow);
//		thermalBoundaryConditions ( flow ); 
//		flow.initializeTemperatures ( tZero );
//		flow.sphericalVertexAreaCalculated=false;
//	}

//	
//        if ( !first && !multithread && (useSolver==0 || fluidBulkModulus>0 || doInterpolate || thermalEngine)) flow.interpolate ( flow.T[!flow.currentTes], flow.tesselation() );
//        if ( waveAction ) flow.applySinusoidalPressure ( flow.tesselation().Triangulation(), sineMagnitude, sineAverage, 30 );
//	else if (boundaryPressure.size()!=0) flow.applyUserDefinedPressure ( flow.tesselation().Triangulation(), boundaryXPos , boundaryPressure);
//        if (normalLubrication || shearLubrication || viscousShear) flow.computeEdgesSurfaces();
//}

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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////....................................Triangulation while maintaining saturation field ............................................//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//void PartialSatClayEngine::reTriangulate()
//{
//      //Governing function to apply triangulation while maintaining saturation distribution.
//      if(debugTPF){std::cerr << endl << "Apply retriangulation";} 
//      initializationTriangulation();
//      readTriangulation();      
//      keepTriangulation = false;
//      initialization();
//      assignWaterVolumesTriangulation();
//      actionMergingAlgorithm();
//      equalizeSaturationOverMergedCells();
//}

//void PartialSatClayEngine::initializationTriangulation()
//{
//  //Resize all relevant functions 
//  
//  //per sphere
//  leftOverVolumePerSphere.resize(scene->bodies->size(),0);
//  untreatedAreaPerSphere.resize(scene->bodies->size(),0);
//  leftOverDVPerSphere.resize(scene->bodies->size(),0);
//  //per tetrahedra
//  finishedUpdating.resize(solver->T[solver->currentTes].cellHandles.size(),0);
//  waterVolume.resize(solver->T[solver->currentTes].cellHandles.size(),0);
//  deltaVoidVolume.resize(solver->T[solver->currentTes].cellHandles.size(),0);
//  tetrahedra.resize(solver->T[solver->currentTes].cellHandles.size());
//  solidFractionSpPerTet.resize(solver->T[solver->currentTes].cellHandles.size());
//  for(unsigned int i =0; i<solver->T[solver->currentTes].cellHandles.size(); i++){
//    tetrahedra[i].resize(4,0);
//    solidFractionSpPerTet[i].resize(4,0);
//  }
//}



//void PartialSatClayEngine::setInitialConditions()
//{
//    if(debugTPF){std::cerr<<endl<<"Set initial condition";}

//    //four possible initial configurations are allowed: primary drainage, primary imbibition, secondary drainage, secondary imbibition
//    RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
//    FiniteCellsIterator cellEnd = tri.finite_cells_end();
//    for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
//        
//    }
//}



//    
//void PartialSatClayEngine::solvePressure()
//{	
//    RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
//    FiniteCellsIterator cellEnd = tri.finite_cells_end();  
//    double oldDT = 0.0;
//    
//    //Define matrix, triplet list, and linear solver
//    tripletList.clear(); // tripletList.resize(T_nnz);
//    Eigen::VectorXd residualsList(numberOfPores);
//    Eigen::VectorXd pressuresList(numberOfPores); //Solve aMatrix * pressuresList = residualsList
//    
//    //define lists
//    if((deformation && remesh) || firstDynTPF){
//	aMatrix.resize(numberOfPores,numberOfPores);
//	saturationList.assign(numberOfPores,0.0);
//	hasInterfaceList.assign(numberOfPores,false);
//	listOfFlux.assign(numberOfPores,0.0);
//	listOfMergedVolume.assign(numberOfPores,0.0);										//NOTE CHANGED AFTER PUSH ON GIT
//    }
//    
//    //reset various lists
//    for(unsigned int i = 0; i < numberOfPores; i++){
//      residualsList[i] = 0.0;
//      pressuresList[i] = 0.0;
//      saturationList[i] = listOfPores[i]->info().saturation;
//      hasInterfaceList[i] = listOfPores[i]->info().hasInterface;
//      listOfFlux[i] = 0.0;
//    } 
//    
//    //Fill matrix
//    for(unsigned int i = 0; i < numberOfPores; i++){
//      
//       //Get diagonal coeff
//       double dsdp2 = 0.0;
//       double coeffA = 0.0, coeffA2 = 0.0;
//       if(hasInterfaceList[i] && !firstDynTPF){
//	    if(listOfPores[i]->info().p() == 0){
//	      std::cout << endl << "Error, pressure = 0 "<< listOfPores[i]->info().p() << listOfPores[i]->info().id;
//	      listOfPores[i]->info().p() = -1.0 * listOfPores[i]->info().thresholdPressure;
//	    }
//	    dsdp2 = dsdp(listOfPores[i],listOfPores[i]->info().p());
//	    coeffA = dsdp2 * ((listOfPores[i]->info().mergedVolume / scene->dt) + (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling)); //Only consider the change in porosity due to particle movement
//	}

//           
//       //fill matrix off-diagonals
//       if(!listOfPores[i]->info().isWResInternal ){	
//	  for(unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++){ 
//	      tripletList.push_back(ETriplet(i,listOfPores[i]->info().poreNeighbors[j],-1.0 * listOfPores[i]->info().listOfkNorm[j]));
//	      coeffA2 += listOfPores[i]->info().listOfkNorm[j];
//	  }
//       }

//	//Set boundary conditions
//	if(listOfPores[i]->info().isWResInternal){	
//	  tripletList.push_back(ETriplet(i,i, 1.0));
//	  residualsList[i] = waterBoundaryPressure;
//	}
//	
//	//Fill matrix diagonal
//	if(!listOfPores[i]->info().isWResInternal ){	
//	  if(hasInterfaceList[i]){residualsList[i] += -1.0*listOfPores[i]->info().saturation * (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling) + coeffA * listOfPores[i]->info().p();}
//	  if(!hasInterfaceList[i] && deformation && listOfPores[i]->info().saturation > listOfPores[i]->info().minSaturation){residualsList[i] += -1.0*(listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling);} 
//	  tripletList.push_back(ETriplet(i,i, coeffA +  coeffA2));
//	}
//    }
//    
//    
//    //Solve Matrix
//    aMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
//    eSolver.analyzePattern(aMatrix);
//    eSolver.factorize(aMatrix);
//    eSolver.compute(aMatrix);
//    
//    
//    //Solve for pressure: FIXME: add check for quality of matrix, if problematic, skip all below. 
//	pressuresList = eSolver.solve(residualsList);
//	
//	//Compute flux
//	double flux = 0.0;
//	double accumulativeDefFlux = 0.0;
//	double waterBefore = 0.0, waterAfter = 0.0;
//	double boundaryFlux = 0.0, lostVolume = 0.0;
//	oldDT = scene->dt;
//	
//	//check water balance
// 	for(unsigned int i = 0; i < numberOfPores; i++){
//	  waterBefore += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume; 
// 	}
//	
//	//compute flux
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  flux = 0.0;
//	    for(unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++){
//		  flux += listOfPores[i]->info().listOfkNorm[j] * (pressuresList[i] - pressuresList[ listOfPores[i]->info().poreNeighbors[j] ]) ;  
//	    }
//	  listOfFlux[i] = flux;
// 	  if(listOfPores[i]->info().saturation > listOfPores[i]->info().minSaturation){accumulativeDefFlux += listOfPores[i]->info().accumulativeDV;}
//	  if(listOfPores[i]->info().isWResInternal){	
//	    boundaryFlux += flux;
//	   }
//	   if(!listOfPores[i]->info().isWResInternal && !hasInterfaceList[i] && std::abs(listOfFlux[i]) > 1e-15 && !deformation){
//	      std::cerr << " | Flux not 0.0" << listOfFlux[i] << " isNWRES:  "<< listOfPores[i]->info().isNWRes << " saturation: "<< listOfPores[i]->info().saturation << " P:"<< listOfPores[i]->info().p() << " isNWef:" << 
//	      listOfPores[i]->info().isNWResDef << "|";
//	      lostVolume += listOfFlux[i] * scene->dt;
//	   }	  	 
//	}

//	
//	double summFluxList = 0.0;
//	double summFluxUnsat = 0.0;
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if(hasInterfaceList[i]){summFluxUnsat += listOfFlux[i];}
//	  summFluxList += listOfFlux[i];
//	}

//	//update saturation
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if(!deformation && hasInterfaceList[i] && listOfFlux[i] != 0.0 && !listOfPores[i]->info().isWResInternal){		
//       	      double ds = -1.0 * scene->dt * (listOfFlux[i]) / (listOfPores[i]->info().mergedVolume +  listOfPores[i]->info().accumulativeDV * scene->dt);
//	      saturationList[i] = ds + (saturationList[i] * listOfPores[i]->info().mergedVolume / (listOfPores[i]->info().mergedVolume +  listOfPores[i]->info().accumulativeDV * scene->dt )); 
//          }
// 	    if(deformation && hasInterfaceList[i] && !listOfPores[i]->info().isWResInternal){
// 		saturationList[i] = saturationList[i] + (pressuresList[i] - listOfPores[i]->info().p())*dsdp(listOfPores[i],listOfPores[i]->info().p()) + (scene->dt / (listOfPores[i]->info().mergedVolume +  (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling) * scene->dt )) * listOfPores[i]->info().accumulativeDVSwelling;
// 	    }
//	  
//	}
//	
//        for(unsigned int i = 0; i < numberOfPores; i++){
//	  waterAfter += saturationList[i] * (listOfPores[i]->info().mergedVolume +   listOfPores[i]->info().accumulativeDV * scene->dt); 							
//        }

//	
//	accumulativeFlux += (summFluxList) * scene->dt;
//	accumulativeDeformationFlux += accumulativeDefFlux * scene->dt;
//	
//	fluxInViaWBC += boundaryFlux * scene->dt;
//	
//	  
// 	if(!deformation && std::abs(boundaryFlux * scene->dt + (waterBefore - waterAfter)) / std::abs(boundaryFlux * scene->dt) > 1e-3 && std::abs(boundaryFlux) > 1e-18){	//FIXME test has to optimized for deforming pore units
// 	  std::cerr << endl << "No volume balance! Flux balance: WBFlux:" << std::abs(boundaryFlux * scene->dt + (waterBefore - waterAfter)) / std::abs(boundaryFlux * scene->dt)  << " " <<boundaryFlux * scene->dt << "Flux: " << summFluxList*scene->dt <<  "deltaVolume: " << waterBefore - waterAfter  << "Flux in IFACE: "<< summFluxUnsat * scene->dt << " lostVolume: " << lostVolume * scene->dt;
//// 	  stopSimulation = true;
// 	}
//	// --------------------------------------find new dt -----------------------------------------------------
//	double dt = 0.0, finalDT = 1e6;
//	int saveID = -1;
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  //Time step for deforming pore units
//	  if(deformation){
//	    dt = -1.0 * listOfPores[i]->info().mergedVolume / (listOfPores[i]->info().accumulativeDV + listOfPores[i]->info().accumulativeDVSwelling);	//Residence time total pore volume
//	    if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;saveID = -1;}
//	    
//	    if(listOfPores[i]->info().accumulativeDVSwelling > 0.0 || listOfPores[i]->info().accumulativeDV > 0.0){ // Residence time during increase in pore size
//	      if(listOfPores[i]->info().accumulativeDVSwelling > listOfPores[i]->info().accumulativeDV){
//		dt = listOfPores[i]->info().mergedVolume * (1.0 - saturationList[i]) / listOfPores[i]->info().accumulativeDVSwelling;
//	      }
//	      if(listOfPores[i]->info().accumulativeDVSwelling <= listOfPores[i]->info().accumulativeDV){
//		dt = listOfPores[i]->info().mergedVolume * (1.0 - saturationList[i]) / listOfPores[i]->info().accumulativeDV;
//	      }
//	      if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;saveID = -2;}
//	    }
//	    if(listOfPores[i]->info().accumulativeDVSwelling < 0.0 || listOfPores[i]->info().accumulativeDV < 0.0){
//		if(listOfPores[i]->info().accumulativeDVSwelling < listOfPores[i]->info().accumulativeDV){
//		  dt = -1.0 * listOfPores[i]->info().mergedVolume * saturationList[i] / listOfPores[i]->info().accumulativeDVSwelling;
//		}
//		if(listOfPores[i]->info().accumulativeDVSwelling >= listOfPores[i]->info().accumulativeDV){
//		  dt = -1.0 * listOfPores[i]->info().mergedVolume * saturationList[i] / listOfPores[i]->info().accumulativeDV;
//		}
//		if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;saveID = -2;}
//	    }
//	  }
//	  
//	  //Time step for dynamic flow
//	  if(hasInterfaceList[i]){		
//	      //thresholdSaturation
//	      if(std::abs(listOfPores[i]->info().thresholdSaturation - saturationList[i]) > truncationPrecision){
//		dt =  -1.0 * (listOfPores[i]->info().thresholdSaturation - saturationList[i]) * listOfPores[i]->info().mergedVolume / listOfFlux[i];
//		if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;/*saveID = 1;*/}
//	      }
//	      //Empty pore
//	      if(std::abs(0.0 - saturationList[i]) > truncationPrecision && listOfFlux[i] > 0.0){ //only for drainage
//		dt =  -1.0 * (0.0 - saturationList[i]) * listOfPores[i]->info().mergedVolume / listOfFlux[i];
//		if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;/*saveID = 2;*/}
//	      }
//	      //Saturated pore
//	      if(std::abs(1.0 - saturationList[i]) > truncationPrecision && listOfFlux[i] < 0.0){ //only for imbibition
//		dt =  -1.0 * (1.0 - saturationList[i]) * listOfPores[i]->info().mergedVolume / listOfFlux[i];
//		if(dt > deltaTimeTruncation && dt < finalDT){finalDT = dt;/*saveID = 3;*/}
//	      }
//	  }
//	}
//	if(finalDT == 1e6){
//            finalDT = deltaTimeTruncation;
//            saveID = 5;
//            if(!firstDynTPF && !remesh){
//                std::cout << endl << "NO dt found!";
//                stopSimulation = true;}            
//        }
//	scene->dt = finalDT * safetyFactorTimeStep;
//	if(debugTPF){std::cerr<<endl<< "Time step: " << finalDT << " Limiting process:" << saveID;}
//	    
//	
//	
//	// --------------------------------------update cappilary pressure (need to correct for linearization of ds/dp)-----------------------------------------------------
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if(hasInterfaceList[i] && !listOfPores[i]->info().isWResInternal && (!deformation || listOfPores[i]->info().saturation != 0.0)){
//	    pressuresList[i] = porePressureFromPcS(listOfPores[i],saturationList[i]);}
//	}
//	

//	// --------------------------------------Find invasion events-----------------------------------------------------
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if(saturationList[i] > 1.0 - truncationPrecision &&  (listOfFlux[i] < 0.0 || deformation) && saturationList[i] != 1.0){
//		if(saturationList[i] > 1.0){
//		 if(saturationList[listOfPores[i]->info().invadedFrom] >= 1.0){
//		   	waterVolumeTruncatedLost +=  (saturationList[i] - 1.0) * listOfPores[i]->info().mergedVolume;
//			saturationList[i] = 1.0;
//		 }
//		 if(saturationList[listOfPores[i]->info().invadedFrom] < 1.0){
//		   saturationList[listOfPores[i]->info().invadedFrom] += (saturationList[i] - 1.0) * listOfPores[i]->info().mergedVolume / listOfPores[listOfPores[i]->info().invadedFrom]->info().mergedVolume;
//		   saturationList[i] = 1.0;
//		   if(saturationList[listOfPores[i]->info().invadedFrom] > 1.0){
//		    waterVolumeTruncatedLost +=  (saturationList[listOfPores[i]->info().invadedFrom] - 1.0) * listOfPores[listOfPores[i]->info().invadedFrom]->info().mergedVolume;
//		    saturationList[listOfPores[i]->info().invadedFrom] = 1.0;
//		   }
//		 }	 
//		 
//		}
// 		saturationList[i] = 1.0;
// 		hasInterfaceList[i] = false;
//		listOfPores[i]->info().isNWResDef = true;
//	  }
//	}
//	
//	//Check for drainage
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if((fractionMinSaturationInvasion == -1 && hasInterfaceList[i] && saturationList[i] <  listOfPores[i]->info().thresholdSaturation) || (fractionMinSaturationInvasion > 0.0 && saturationList[i] < fractionMinSaturationInvasion)){
//	    for(unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++){
//	      if(airBoundaryPressure - pressuresList[listOfPores[i]->info().poreNeighbors[j]] > listOfPores[i]->info().listOfEntryPressure[j] &&
//		!hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] && 
//		!listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isWResInternal && 
//		saturationList[listOfPores[i]->info().poreNeighbors[j]] > listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().minSaturation &&
//		saturationList[listOfPores[i]->info().poreNeighbors[j]] <= 1.0		
//	      ){
//		  hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] = true;
//		  saturationList[listOfPores[i]->info().poreNeighbors[j]] = 1.0 - truncationPrecision;
//		  pressuresList[listOfPores[i]->info().poreNeighbors[j]] = porePressureFromPcS(listOfPores[i], 1.0 - truncationPrecision);
//		  listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().invadedFrom = i;

// 	      }
//	    }
//	  }
//	}
//	    
//	
//	
//	//truncate saturation
//	for(unsigned int i = 0; i < numberOfPores; i++){
//	  if((saturationList[i] < truncationPrecision || saturationList[i] <= listOfPores[i]->info().minSaturation)){		
//	    waterVolumeTruncatedLost -=  (listOfPores[i]->info().minSaturation - saturationList[i]) * listOfPores[i]->info().mergedVolume;
//	    saturationList[i] = listOfPores[i]->info().minSaturation; 
//// 	    hasInterfaceList[i] = false; // NOTE: in case of deactivation of empty cell, set hasInterfaceList[i] to false
//	    pressuresList[i] =  porePressureFromPcS(listOfPores[i], listOfPores[i]->info().minSaturation); //waterBoundaryPressure;				
//	    listOfPores[i]->info().isNWRes = true;
//	    for(unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++){
//		if(!hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] && !listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isNWRes && listOfFlux[listOfPores[i]->info().poreNeighbors[j]] >= 0.0 && !listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isWResInternal){
//		  hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] = true;
//		  saturationList[listOfPores[i]->info().poreNeighbors[j]] = 1.0 - truncationPrecision;
//		  pressuresList[listOfPores[i]->info().poreNeighbors[j]] = porePressureFromPcS(listOfPores[listOfPores[i]->info().poreNeighbors[j]], 1.0 - truncationPrecision);
//		}
//	    }
//	  }
//	  
//	  if(listOfPores[i]->info().isNWResDef && saturationList[i] < listOfPores[i]->info().thresholdSaturation){
//	    listOfPores[i]->info().isNWResDef = false;
//	  }
//	  
//	}
// 



//    
//    if(deformation){
//	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++){
//	  cell->info().poreBodyVolume += cell->info().dv() * oldDT;			
//	}
//	//copyPoreDataToCells();  //NOTE: For two-way coupling this function should be activated, but it is a bit costly for computations.
//    }
//    for(unsigned int i = 0; i < numberOfPores; i++){
//		if(deformation){
//		  listOfPores[i]->info().mergedVolume  += listOfPores[i]->info().accumulativeDV * oldDT;
//		  listOfPores[i]->info().poreBodyRadius = getChi(listOfPores[i]->info().numberFacets)*std::pow(listOfPores[i]->info().mergedVolume,(1./3.));			
//		}
//		if(saturationList[i]  > 1.0){
//		 std::cerr << endl << "Error!, saturation larger than 1? "; 
//		 saturationList[i] = 1.0;								//NOTE ADDED AFTER TRUNK UPDATE should be 0.0?
//// 		 stopSimulation = true;
//		}
//		listOfPores[i]->info().saturation = saturationList[i];
//		listOfPores[i]->info().p() = pressuresList[i];
//		listOfPores[i]->info().hasInterface = bool(hasInterfaceList[i]);
//		listOfPores[i]->info().flux = listOfFlux[i];
//		listOfPores[i]->info().dv() = 0.0;							//NOTE ADDED AFTER TRUNK UPDATE
//		listOfPores[i]->info().accumulativeDV = 0.0;
//    }
//}


//        
//void PartialSatClayEngine::actionTPF()
//{
//  iterationTPF += 1;
//  if(firstDynTPF){
//      std::cout << endl << "Welcome to the two-phase flow Engine" << endl << "by T.Sweijen, B.Chareyre and S.M.Hassanizadeh" << endl << "For contact: T.Sweijen@uu.nl";
//      solver->computePermeability();
//      scene->time = 0.0;
//      initialization();
//      actionMergingAlgorithm();
//      calculateResidualSaturation();
//      setInitialConditions();
//      setBoundaryConditions();
//      verifyCompatibilityBC();
//      setPoreNetwork();
//      scene->dt = 1e-20;
//      setListOfPores();
//      solvePressure();
//      getQuantities();
//      firstDynTPF = false;
//  }
//  if(!firstDynTPF && !stopSimulation){
////    bool remesh = false;
//   //Time steps + deformation, but no remeshing
//   scene->time = scene->time + scene->dt;
//   if(deformation && !remesh){
//     updateDeformationFluxTPF();
//     if(int(float(iterationTPF)/10.0) == float(iterationTPF)/10.0){
//	updatePoreUnitProperties();
//   }
//  }
//   //Update pore throat radii etc. 
//   
//   if(deformation && remesh){
//     reTriangulate(); //retriangulation + merging
//     calculateResidualSaturation();
//     transferConditions(); //get saturation, hasInterface from previous network
//     setBoundaryConditions();
//     setPoreNetwork();
//   }
//   setListOfPores();
//   if(solvePressureSwitch){solvePressure();}
//   if(deformation){if(int(float(iterationTPF)/50.0) == float(iterationTPF)/50.0){getQuantities();}}	//FIXME update of quantities has to be made more appropiate

//   
////    getQuantities();//NOTE FIX
//   
//    if(!deformation ){
//      if(!getQuantitiesUpdateCont){
//	if(int(float(iterationTPF)/100.0) == float(iterationTPF)/100.0){getQuantities();}
//      }
//      if(getQuantitiesUpdateCont){
//	getQuantities();
//      }
//    }

//   
//   if(remesh){remesh = false;} //Remesh bool is also used in solvePressure();

//  }
//}


#endif //PartialSat


