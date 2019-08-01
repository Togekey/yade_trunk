// *  Deepak kn : deepak.kunhappan@3sr-grenoble.fr; deepak.kn1990@gmail.com

#ifdef YADE_MPI

#include <mpi.h>
#include "FoamCoupling.hpp"

#include<pkg/common/Facet.hpp>
#include<pkg/common/Box.hpp>
#include<pkg/common/Sphere.hpp>
#include<pkg/common/Grid.hpp>

#include <iostream>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(FoamCoupling);
YADE_PLUGIN((FoamCoupling));
YADE_PLUGIN((FluidDomainBbox));
CREATE_LOGGER(FluidDomainBbox);
YADE_PLUGIN((Bo1_FluidDomainBbox_Aabb));


void Bo1_FluidDomainBbox_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* b){
	
	FluidDomainBbox* domain = static_cast<FluidDomainBbox*>(cm.get());
	if (!bv){bv = shared_ptr<Bound> (new Aabb); }
	Aabb* aabb = static_cast<Aabb*>(bv.get()); 
	aabb-> min = domain -> minBound; 
	aabb-> max = domain -> maxBound; 
}


void FoamCoupling::getRank() {

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

}

void FoamCoupling::setNumParticles(int np){

      getRank(); 
      numParticles = np;
      castNumParticle(numParticles); 
      initDone = true; 
  
}

void FoamCoupling::setIdList(const std::vector<int>& alist) {
  bodyList.clear(); bodyList.resize(alist.size()); 
  for (unsigned int i=0; i != bodyList.size(); ++i){
    bodyList[i] = alist[i];
  }

}


void FoamCoupling::insertBodyId(int bId){
	bodyList.push_back(bId);
  
}

bool FoamCoupling::eraseId(int bId){
	auto it = std::find(bodyList.begin(), bodyList.end(), bId);
	if (it != bodyList.end()){bodyList.erase(it); return true; }
	else {
		LOG_ERROR("Id not found in list of ids in coupling"); 
		return false; 
	}
}


int FoamCoupling::getNumBodies(){
	return bodyList.size(); 
}

std::vector<int> FoamCoupling::getIdList(){
	return bodyList; 
}

void FoamCoupling::castParticle() {
  
	int sz = bodyList.size(); 
	MPI_Bcast(&sz, 1, MPI_INT, rank, MPI_COMM_WORLD);
	procList.resize(sz); hydroForce.resize(sz*6); 
	particleData.resize(10*sz); 
	std::fill(procList.begin(), procList.end(), -1); 
	std::fill(hydroForce.begin(), hydroForce.end(), 1e-50); 

#ifdef YADE_OPENMP
#pragma omp parallel  for collapse (1)
#endif

for (unsigned int i=0; i <  bodyList.size(); ++i)
  {
    const Body* b = (*scene -> bodies)[bodyList[i]].get();
    if ( scene-> isPeriodic){
      Vector3r pos = scene->cell->wrapPt( b->state->pos);
      particleData[i*10] = pos[0];
      particleData[i*10+1] = pos[1];
      particleData[i*10+2] = pos[2];
    } else {

    particleData[i*10] = b->state->pos[0];
    particleData[i*10+1] = b->state->pos[1];
    particleData[i*10+2] = b->state->pos[2];
  }
    particleData[i*10+3] = b->state->vel[0];
    particleData[i*10+4] = b->state->vel[1];
    particleData[i*10+5] = b->state->vel[2];
    particleData[i*10+6] = b->state->angVel[0];
    particleData[i*10+7] = b->state->angVel[1];
    particleData[i*10+8] = b->state->angVel[2];
    shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
    particleData[i*10+9] = s->radius;
  }
  MPI_Bcast(&particleData.front(), particleData.size(), MPI_DOUBLE, rank, MPI_COMM_WORLD);
  // clear array after bcast
  particleData.clear(); 

}

void FoamCoupling::castNumParticle(int value) {
	MPI_Bcast(&value, 1, MPI_INT, rank, MPI_COMM_WORLD);
}


void FoamCoupling::castTerminate() {
  int value = 10; 
  MPI_Bcast(&value, 1, MPI_INT, rank, MPI_COMM_WORLD);

}

void FoamCoupling::updateProcList()
{
  for (unsigned int i=0; i != bodyList.size(); ++i)
  {
     int dummy_val = -5;
     MPI_Allreduce(&dummy_val,&procList[i],1,MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     if (procList[i] < 0 )  std::cout << "Particle not found in FOAM " << std::endl;
   }
   
   
}

void FoamCoupling::recvHydroForce() {
  for (unsigned int i=0; i!= procList.size(); ++i) {
    int recvFrom = procList[i];
    for (unsigned int j=0; j != 6; ++j) {
     MPI_Recv(&hydroForce[6*i+j],1,MPI_DOUBLE,recvFrom,sendTag,MPI_COMM_WORLD,&status); 
    }
     
}}

void FoamCoupling::setHydroForce() {

  // clear hydroforce before summation
  
 #ifdef YADE_OPENMP
 #pragma omp parallel for collapse(1)
 #endif
    for (unsigned int i=0; i < bodyList.size(); ++i) {
       const Vector3r& fx=Vector3r(hydroForce[6*i], hydroForce[6*i+1], hydroForce[6*i+2]);
       const Vector3r& tx=Vector3r(hydroForce[6*i+3], hydroForce[6*i+4], hydroForce[6*i+5]);
       scene->forces.addForce(bodyList[i], fx);
       scene->forces.addTorque(bodyList[i], tx);
  }

}

void FoamCoupling::sumHydroForce() {
  // clear the vector
  std::fill(hydroForce.begin(), hydroForce.end(), 0.0);
  Real dummy_val = 0.0;
  for (unsigned int i=0; i != bodyList.size(); ++i) {
   for (unsigned int j=0; j != 6; ++j){
     MPI_Allreduce(&dummy_val ,&hydroForce[6*i+j],1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
  }
}


void FoamCoupling::resetProcList() {
     procList.clear(); 
}


void FoamCoupling::action() {

	if ( exchangeData()) {
		runCoupling();
		exchangeDeltaT();
	}
	setHydroForce();
}

bool FoamCoupling::exchangeData(){

	return scene->iter%dataExchangeInterval==0;

}

void FoamCoupling::exchangeDeltaT() {

  // Recv foamdt  first and broadcast;
  MPI_Recv(&foamDeltaT,1,MPI_DOUBLE,1,sendTag,MPI_COMM_WORLD,&status);
  //bcast yadedt to others.
  Real  yadeDt = scene-> dt;
  MPI_Bcast(&yadeDt,1,MPI_DOUBLE, rank, MPI_COMM_WORLD);
  // calculate the interval . TODO: to include hydrodynamic time scale if inertial in openfoam
  // here -> hDeltaT = getViscousTimeScale();
  dataExchangeInterval = (long int) ((yadeDt < foamDeltaT) ? foamDeltaT/yadeDt : 1);

}

Real FoamCoupling::getViscousTimeScale() {

//  Real hDeltaT = 0.0;
//  Real dummy = 1e9;
//
//  MPI_Allreduce(&dummy, &hDeltaT, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//
	return 0;
}

void FoamCoupling::runCoupling() {
	castParticle();
	updateProcList();
	if (isGaussianInterp){
		sumHydroForce();
	} else {
		recvHydroForce();
	}
}


void FoamCoupling::getFluidDomainBbox() {
	
	
	int localCommSize, localRank, worldCommSize, worldRank; 
	const shared_ptr<Body>& thisSubdomainBody = (*scene->bodies)[scene->thisSubdomainId]; 
	const shared_ptr<Subdomain>& thisSd = YADE_PTR_CAST<Subdomain>(thisSubdomainBody->shape); 
	
	
	//get local comm size
	MPI_Comm_size(thisSd->selfComm(), &localCommSize); 
	MPI_Comm_rank(thisSd->selfComm(), &localRank); 
	
	//world comm size and world rank 
	
	MPI_Comm_size(MPI_COMM_WORLD, &worldCommSize); 
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank); 
	
	commSzdff = abs(localCommSize-worldCommSize); 
	//std::vector<int> fluidRanks(szdff, -1); 
	
	std::vector<std::vector<double> > minMaxBuff; 
	
	//alloc memory 
	for (int i=0; i != commSzdff; ++i){
		std::vector<double> buff(6, 1e-50); 
		minMaxBuff.push_back(std::move(buff)); 
	  
	}

	//recv thegird minmax from fluid solver. 
	for (int rnk=0; rnk != commSzdff; ++rnk){
		std::vector<double>& buff = minMaxBuff[rnk]; 
		MPI_Bcast(&buff.front(), 6,MPI_DOUBLE, rnk+localCommSize, MPI_COMM_WORLD); 
		
	}
	
	//create fluidDomainBbox bodies and get their ids. 
	for (int fd = 0; fd != commSzdff; ++fd){
		shared_ptr<Body>  flBody(shared_ptr<Body> (new Body()));
		shared_ptr<FluidDomainBbox> flBodyshape(shared_ptr<FluidDomainBbox>  (new FluidDomainBbox())); 
		std::vector<double>& buff = minMaxBuff[fd];
		flBodyshape->setMinMax(buff); 
		flBodyshape->domainRank = localCommSize+fd; 
		flBodyshape->hasIntersection = false; 
		flBody->shape = flBodyshape;  // casting?!  
		fluidDomains.push_back(scene->bodies->insert(flBody)); 
		
	}
	
  
}

void FoamCoupling::buildSharedIds() {
	const shared_ptr<Body>& subdBody = (*scene->bodies)[scene->thisSubdomainId]; 
	const shared_ptr<Subdomain>& subD = YADE_PTR_CAST<Subdomain>(subdBody->shape); 
	for (int bodyId = 0; bodyId != static_cast<int>(subD->ids.size()); ++bodyId){
		std::vector<Body::id_t> fluidIds; 
		const shared_ptr<Body>& testBody = (*scene->bodies)[subD->ids[bodyId]]; 
		for (const auto& bIntrs : testBody->intrs){
			const shared_ptr<Interaction>& intr = bIntrs.second; 
			Body::id_t otherId; 
			if (testBody->id == intr->getId1()){otherId = intr->getId2();  } else { otherId = intr->getId2(); }
			if (ifFluidDomain(otherId)){ fluidIds.push_back(otherId); }
		}
		if (fluidIds.size()){sharedIds.push_back(std::make_pair(subD->ids[bodyId], fluidIds)); } 
	}
}

unsigned  FoamCoupling::ifSharedId(const Body::id_t& testId){
	
	if (! sharedIds.size()) {return false; } 
	bool inIds; unsigned res = 0; 
	for (unsigned indx =0; indx != sharedIds.size(); ++indx ){
		if (testId == sharedIds[indx].first) {
			res = testId; 
		}
	}
	return res; 
	
}


// bool FoamCoupling::


bool FoamCoupling::ifFluidDomain(const Body::id_t&  testId ){ 
	const auto& iter = std::find(fluidDomains.begin(), fluidDomains.end(), testId);  
	if (iter != fluidDomains.end()){return true; } else {return false;}
  
}

void FoamCoupling::findIntersections(){
	
	/*find bodies intersecting with the fluid domain bbox, get their ids and ranks of the intersecting fluid processes.  */
		
	for (unsigned f = 0; f != fluidDomains.size(); ++f){
		shared_ptr<Body>& fdomain = (*scene->bodies)[fluidDomains[f]]; 
		if (fdomain){
			for (const auto& fIntrs : fdomain->intrs){
				Body::id_t otherId; 
				const shared_ptr<Interaction>& intr = fIntrs.second;
				if (fdomain->id == intr->getId1()){otherId = intr->getId2(); } else { otherId = intr->getId1(); }
				const shared_ptr<Body>& testBody = (*scene->bodies)[otherId]; 
				if (testBody) {
					if (testBody->subdomain==scene->subdomain){
						if (!ifDomainBodies(testBody)){
							const shared_ptr<FluidDomainBbox>& flBox = YADE_PTR_CAST<FluidDomainBbox>(fdomain->shape); 
							flBox->bIds.push_back(testBody->id); 
							flBox->hasIntersection = true; 
						}
					}
				} 
			}
		
	}
	}
}

bool FoamCoupling::ifDomainBodies(const shared_ptr<Body>& b) {
	// check if body is subdomain, wall, facet, or other fluidDomainBbox 
	
	shared_ptr<Box> boxShape = YADE_PTR_DYN_CAST<Box> (b->shape); 
	shared_ptr<FluidDomainBbox> fluidShape = YADE_PTR_DYN_CAST<FluidDomainBbox> (b->shape); 
	shared_ptr<Facet> facetShape = YADE_PTR_DYN_CAST<Facet>(b->shape); 
	
	if (b->getIsSubdomain()){return true; }
	else if (boxShape) {return true; }
	else if (facetShape) {return true; }
	else {return false; }
	 
}


void FoamCoupling::sendIntersectionToFluidProcs(){
	
	// notify the fluid procs about intersection based on number of intersecting bodies. 
	// vector of sendRecvRanks, with each vector element containing the number of bodies, if no bodies
	sendRecvRanks.resize(fluidDomains.size()); 
	
	for (unsigned f=0;  f != fluidDomains.size(); ++f){
		const shared_ptr<Body>& fdomain = (*scene->bodies)[fluidDomains[f]]; 
		if (fdomain){
			const shared_ptr<FluidDomainBbox>& fluidBox = YADE_PTR_CAST<FluidDomainBbox>(fdomain->shape); 
			if (fluidBox->hasIntersection){
				sendRecvRanks[f] = fluidBox->bIds.size(); // how about setting up of ids& data here? 
				
			} else {sendRecvRanks[f] = -1; }
		} 
		else {sendRecvRanks[f] = -1; }
	}
	// 
	int buffSz = fluidDomains.size(); 
	
	//MPI_Send or bcast.. 
	for (int rnk = 0; rnk != commSzdff; ++ rnk){
		MPI_Send(&sendRecvRanks.front(), buffSz, MPI_INT, rnk+commSzdff, 2500, MPI_COMM_WORLD); 
		
	}
	
	
	
}

void FoamCoupling::sendBodyData(){
	/* send the particle data to the associated fluid procs. prtData -> pos, vel, angvel, raidus (for sphere), if fiber -> ori  */ 
  
	bool isPeriodic = scene->isPeriodic; 
	for (int f = 0; f != static_cast<int>(fluidDomains.size()); ++f){
		const shared_ptr<Body>& flbody = (*scene->bodies)[fluidDomains[f]];
		if (flbody){
		const shared_ptr<FluidDomainBbox> flbox = YADE_PTR_CAST<FluidDomainBbox>(flbody->shape); 
			if (flbox->hasIntersection){
				std::vector<double> prtData(10*flbox->bIds.size(), 1e-50);
				for (unsigned int i=0; i != flbox->bIds.size(); ++i){
					const shared_ptr<Body>& b = (*scene->bodies)[flbox->bIds[i]]; 
					if (isPeriodic){
						const Vector3r& pos = scene->cell->wrapPt(b->state->pos); 
						prtData[10*i] = pos[0]; 
						prtData[10*i+1] = pos[1]; 
						prtData[10*i+2] = pos[2]; 
						
					} else {
						prtData[10*i] = b->state->pos[0]; 
						prtData[10*i+1] = b->state->pos[1];
						prtData[10*i+2] = b->state->pos[2]; 
					}
					prtData[10*i+3] = b->state->vel[0]; 
					prtData[10*i+4] = b->state->vel[1]; 
					prtData[10*i+5] = b->state->vel[2]; 
					prtData[10*i+6] = b->state->angVel[0]; 
					prtData[10*i+7] = b->state->angVel[1]; 
					prtData[10*i+8] = b->state->angVel[2]; 
					
					const shared_ptr<Sphere>& sph = YADE_PTR_CAST<Sphere>(b->shape); 
					prtData[10*i+9] = sph->radius; 
				}
				
				int sz = prtData.size();  
				MPI_Send(&prtData.front(),sz, MPI_DOUBLE, flbox->domainRank, 1000, MPI_COMM_WORLD ); 
			}
		}
	}
}


void FoamCoupling::verifyParticleDetection() {
	//TODO: from here
// 	std::map<int, std::vector<int> > verifyTracking; 
// 	for (unsigned int  f=0; f != fluidDomains.size(); ++f) {
// 		const shared_ptr<Body>& flbdy  = (*scene->bodies)[fluidDomains[f]]; 
// 		if (flbdy){
// 			const shared_ptr<FluidDomainBbox>& flbox = YADE_PTR_CAST<FluidDomainBbox>(flbdy->shape);
// 			std::vector<int> vt(flbox->bIds.size(), -1); 
// 			verifyTracking.insert(std::make_pair(flbox->domainRank, std::move(vt))); 
// 		}
// 	}
// 	
// 	for (auto& it : verifyTracking){
// 		std::vector<int>& vt = it.second; 
// 		int rnk = it.first; 
// 		MPI_Status status; 
// 		int buffSz = vt.size(); 
// 		MPI_Recv(&vt.front(), buffSz, MPI_INT, rnk , 1002, MPI_COMM_WORLD, &status); 
// 		
// 	}
// 	
// 
// 	for (const auto& it : verifyTracking){
// 		const std::vector<int>& testV = it.second; 
// 		int indx = abs((it.first)-commSzdff); 
// // 		
// 		
// 	}
// 		
// 		
	
	
	
}


bool FoamCoupling::checkSharedDomains(const int& indx) {
	const std::vector<int>& fDomains = sharedIds[indx].second; 
	return true; 
  
}


void FoamCoupling::getParticleForce(){
	return ; 
}



void FoamCoupling::killMPI() { 
	castTerminate(); 
	MPI_Finalize();

}


#endif
}
