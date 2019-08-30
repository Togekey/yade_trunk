#include <core/ForceContainer.hpp>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<core/BodyContainer.hpp>
#include<core/Scene.hpp>

CREATE_LOGGER(ForceContainer);

void ForceContainer::ensureSynced() {
  if(!synced) throw runtime_error("ForceContainer not thread-synchronized; call sync() first!");
}

void ForceContainer::addForceUnsynced(Body::id_t id, const Vector3r& f) {
  assert ((size_t)id<size);
  _force[id]+=f;
}

void ForceContainer::addTorqueUnsynced(Body::id_t id, const Vector3r& m) {
  assert ((size_t)id<size);
  _torque[id]+=m;
}

void ForceContainer::resizePerm(size_t newSize) {
  if (newSize<_permForce.size()) LOG_WARN("permForce may have been assigned to an id larger than maxId, and will be ignored in that case");
  if (newSize>_permForce.size()){
    _permForce.reserve(size_t(1.3*newSize));
    _permTorque.reserve(size_t(1.3*newSize));
    _permForce.resize(newSize,Vector3r::Zero());
    _permTorque.resize(newSize,Vector3r::Zero());
    syncedSizes=false;}
}

#ifdef YADE_OPENMP
#include <omp.h>
void ForceContainer::ensureSize(Body::id_t id, int threadN) {
// 	assert(id < size);
// 	if (not (id < Body::id_t(size))) LOG_ERROR("ensureSize failed "<<threadN);
// 	std::cerr<<"ensureZie(1)"<< Omega::instance().getScene()->subdomain<<" "<<Omega::instance().getScene()->iter <<" "<< id<<" "<<threadN<<std::endl;
//   assert(nThreads>omp_get_thread_num());
//   const Body::id_t idMaxTmp = max(id, _maxId[threadN]);
	_maxId[threadN] = std::max( _maxId[threadN],id);
// //   std::cerr<<"ensureZie(2) "<< Omega::instance().getScene()->subdomain<<" "<<Omega::instance().getScene()->iter <<" "<< threadN<<std::endl;
  if (threadN<0) resizePerm(id+1);
  else if (not (sizeOfThreads[threadN]>unsigned(_maxId[threadN]))) resize(_maxId[threadN]+1,threadN);
  
//   std::cerr<<"ensureZie(3) "<< Omega::instance().getScene()->subdomain<<" "<<Omega::instance().getScene()->iter <<" "<< threadN<<std::endl;
}

ForceContainer::ForceContainer() {
  nThreads=omp_get_max_threads();
  size=0;
  for(int i=0; i<nThreads; i++){
    _forceData.push_back(vvector());
    _torqueData.push_back(vvector());
    _moveData.push_back(vvector());
    _rotData.push_back(vvector());
    sizeOfThreads.push_back(0);
    _maxId.push_back(0);
  }
}

const Vector3r& ForceContainer::getForce(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_force[id]:_zero;
}

void ForceContainer::addForce(Body::id_t id, const Vector3r& f){
//   ensureSize(id,omp_get_thread_num());
   synced=false;
  _forceData[omp_get_thread_num()][id]+=f;
}

const Vector3r& ForceContainer::getTorque(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_torque[id]:_zero;
}

void ForceContainer::addTorque(Body::id_t id, const Vector3r& t) {
//   ensureSize(id,omp_get_thread_num());
   synced=false;
  _torqueData[omp_get_thread_num()][id]+=t;
}

const Vector3r& ForceContainer::getMove(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_move[id]:_zero;
}

void ForceContainer::addMove(Body::id_t id, const Vector3r& m) {
  ensureSize(id,omp_get_thread_num());
  synced=false;
  moveRotUsed=true;
  _moveData[omp_get_thread_num()][id]+=m;
}

const Vector3r& ForceContainer::getRot(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_rot[id]:_zero;
}

void ForceContainer::addRot(Body::id_t id, const Vector3r& r) {
  ensureSize(id,omp_get_thread_num());
  synced=false;
  moveRotUsed=true;
  _rotData[omp_get_thread_num()][id]+=r;
}

void ForceContainer::addMaxId(Body::id_t id) {
	_maxId[omp_get_thread_num()]=std::max(id,_maxId[omp_get_thread_num()]);
}

void ForceContainer::setPermForce(Body::id_t id, const Vector3r& f) {
  ensureSize(id,-1);
  synced=false;
  _permForce[id]=f;
  permForceUsed=true;
}

void ForceContainer::setPermTorque(Body::id_t id, const Vector3r& t) {
  ensureSize(id,-1);
  synced=false;
  _permTorque[id]=t;
  permForceUsed=true;
}

const Vector3r& ForceContainer::getPermForce(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_permForce[id]:_zero;
}

const Vector3r& ForceContainer::getPermTorque(Body::id_t id) {
  ensureSynced();
  return ((size_t)id<size)?_permTorque[id]:_zero;
}

const Vector3r& ForceContainer::getForceUnsynced(Body::id_t id) {
  assert ((size_t)id<size);
  return _force[id];
}

const Vector3r& ForceContainer::getTorqueUnsynced(Body::id_t id) {
  assert ((size_t)id<size);
  return _torque[id];
}

const Vector3r ForceContainer::getForceSingle(Body::id_t id) {
  Vector3r ret(Vector3r::Zero());
  for(int t=0; t<nThreads; t++) {
    ret+=((size_t)id<sizeOfThreads[t])?_forceData [t][id]:_zero;
  }
  if (permForceUsed) ret+=_permForce[id];
  return ret;
}

const Vector3r ForceContainer::getTorqueSingle(Body::id_t id) {
  Vector3r ret(Vector3r::Zero());
  for(int t=0; t<nThreads; t++) {
    ret+=((size_t)id<sizeOfThreads[t])?_torqueData[t][id]:_zero;
  }
  if (permForceUsed) ret+=_permTorque[id];
  return ret;
}

const Vector3r ForceContainer::getMoveSingle(Body::id_t id) {
  Vector3r ret(Vector3r::Zero());
  for(int t=0; t<nThreads; t++) {
    ret+=((size_t)id<sizeOfThreads[t])?_moveData[t][id]:_zero;
  }
  return ret;
}

const Vector3r ForceContainer::getRotSingle(Body::id_t id) {
  Vector3r ret(Vector3r::Zero());
  for(int t=0; t<nThreads; t++) {
    ret+=((size_t)id<sizeOfThreads[t])?_rotData[t][id]:_zero;
  }
  return ret;
}

void ForceContainer::sync(){
  if(synced) return;
  boost::mutex::scoped_lock lock(globalMutex);
  if(synced) return; // if synced meanwhile

  syncSizesOfContainers();
  
#ifdef YADE_MPI
  Omega::instance().getScene()->bodies->updateSubdomainLists();
  const vector<Body::id_t>& boundedSubDBodies = Omega::instance().getScene()->bodies->boundedSubDBodies;
  const unsigned long len=(long)boundedSubDBodies.size();  
  #pragma omp parallel for schedule(static)
  for(unsigned long k=0; k<len; k++){
	  const Body::id_t& id=boundedSubDBodies[k];
#else
  for(long id=0; id<(long)size; id++){
#endif
    Vector3r sumF(Vector3r::Zero()), sumT(Vector3r::Zero());
    for(int thread=0; thread<nThreads; thread++){ sumF+=_forceData[thread][id]; sumT+=_torqueData[thread][id];
	    _forceData[thread][id]=Vector3r::Zero(); _torqueData[thread][id]=Vector3r::Zero(); }  //reset here so we don't have to do it later
    _force[id]+=sumF; _torque[id]+=sumT;
    if (permForceUsed) {_force[id]+=_permForce[id]; _torque[id]+=_permTorque[id];}
  }
  if(moveRotUsed){
    for(long id=0; id<(long)size; id++){
      Vector3r sumM(Vector3r::Zero()), sumR(Vector3r::Zero());
      for(int thread=0; thread<nThreads; thread++){ sumM+=_moveData[thread][id]; sumR+=_rotData[thread][id];}
      _move[id]=sumM; _rot[id]=sumR;
    }
  }
  synced=true; syncCount++;
}

#ifndef YADE_MPI

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
// this is to remove warning about manipulating raw memory
#pragma GCC diagnostic ignored "-Wclass-memaccess"
void ForceContainer::reset(long iter, bool resetAll) {
  syncSizesOfContainers();
  for(int thread=0; thread<nThreads; thread++){
    memset(&_forceData [thread][0],0,sizeof(Vector3r)*sizeOfThreads[0]);
    memset(&_torqueData[thread][0],0,sizeof(Vector3r)*sizeOfThreads[0]);
    if(moveRotUsed){
      memset(&_moveData  [thread][0],0,sizeof(Vector3r)*sizeOfThreads[0]);
      memset(&_rotData   [thread][0],0,sizeof(Vector3r)*sizeOfThreads[0]);
    }
  }
  memset(&_force [0], 0,sizeof(Vector3r)*size);
  memset(&_torque[0], 0,sizeof(Vector3r)*size);
  if(moveRotUsed){
    memset(&_move  [0], 0,sizeof(Vector3r)*size);
    memset(&_rot   [0], 0,sizeof(Vector3r)*size);
  }
  if (resetAll and permForceUsed){
    memset(&_permForce [0], 0,sizeof(Vector3r)*size);
    memset(&_permTorque[0], 0,sizeof(Vector3r)*size);
    permForceUsed = false;
  }
  if (!permForceUsed) synced=true; else synced=false;
  moveRotUsed=false;
  lastReset=iter;
}

#else

int w=0;
void ForceContainer::reset(long iter, bool resetAll) {
	syncSizesOfContainers();
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	scene->bodies->updateSubdomainLists();
	const vector<Body::id_t>& sdIds = scene->bodies->boundedSubDBodies;
	
	size_t currSize=sdIds.size();
	
	#pragma omp parallel for schedule(static)
	for(int thread=0; thread<nThreads; thread++){
		// no need to reset threads since they are set to zero in sync() already
// 		for (unsigned long k=0;k<currSize;k++) _forceData[thread][sdIds[k]]= Vector3r::Zero();
// 		if (scene->subdomain==w) std::cerr<<"("<<scene->subdomain<<")AB: "<<thread<<"  "<<sizeOfThreads[thread]<<" "<<sdIds.size()<<std::endl;
// 		for (unsigned long k=0;k<currSize;k++) _torqueData[thread][sdIds[k]]=Vector3r::Zero();
// 		if (scene->subdomain==w) std::cerr<<"("<<scene->subdomain<<")B: "<<thread<<" "<<sizeOfThreads[thread]<<" "<<sdIds.size()<<std::endl;
		if(moveRotUsed){
			for (unsigned long k=0;k<currSize;k++) _moveData[thread][sdIds[k]]=Vector3r::Zero();
			for (unsigned long k=0;k<currSize;k++) _rotData[thread][sdIds[k]]=Vector3r::Zero();
		}
	}
	#pragma omp parallel for schedule(static)
	for (unsigned long k=0;k<currSize;k++) _force[sdIds[k]]=Vector3r::Zero();
	#pragma omp parallel for schedule(static)
	for (unsigned long k=0;k<currSize;k++) _torque[sdIds[k]]=Vector3r::Zero();
	if(moveRotUsed){
		for (unsigned long k=0;k<size;k++) _move[k]=Vector3r::Zero();
		for (unsigned long k=0;k<size;k++) _force[k]=Vector3r::Zero();
	}
	if (resetAll){
		for (unsigned long k=0;k<currSize;k++) _permForce[sdIds[k]]=Vector3r::Zero();
		for (unsigned long k=0;k<currSize;k++) _permTorque[sdIds[k]]=Vector3r::Zero();
		permForceUsed = false;
	}
	if (!permForceUsed) synced=true; else synced=false;
	moveRotUsed=false;
	lastReset=iter;
}

#endif
#pragma GCC diagnostic pop

void ForceContainer::resize(size_t newSize, int threadN) {
  if (sizeOfThreads[threadN]>=newSize) return;
  LOG_DEBUG("Resize ForceContainer from the size "<<size<<" to the size "<<newSize);
  _forceData[threadN].reserve(size_t(newSize*1.3));
  _torqueData[threadN].reserve(size_t(newSize*1.3));
  _forceData[threadN].resize(newSize,Vector3r::Zero());
  _torqueData[threadN].resize(newSize,Vector3r::Zero());
  if (moveRotUsed) {
    _moveData[threadN].reserve(size_t(newSize*1.3));
    _rotData[threadN].reserve(size_t(newSize*1.3));
    _moveData[threadN].resize(newSize,Vector3r::Zero());
    _rotData[threadN].resize(newSize,Vector3r::Zero());}
  sizeOfThreads[threadN] = newSize;
  _maxId[threadN]=newSize-1;
  syncedSizes=false;
}

int ForceContainer::getNumAllocatedThreads() const {return nThreads;}
bool ForceContainer::getMoveRotUsed() const {return moveRotUsed;}
bool ForceContainer::getPermForceUsed() const {return permForceUsed;}

void ForceContainer::syncSizesOfContainers() {
  //check whether all containers have equal length, and if not resize it  
  size_t maxThreadSize = 0;
  for(int i=0; i<nThreads; i++) maxThreadSize = std::max(maxThreadSize,size_t(_maxId[i]+1));
  if (maxThreadSize>size) syncedSizes=false;
  if (syncedSizes) return;
  size_t newSize = std::max(size,maxThreadSize);  
  for(int i=0; i<nThreads; i++) resize(newSize,i);

  if (newSize>size){
	  _force.reserve(size_t(newSize*1.3));
	  _torque.reserve(size_t(newSize*1.3));
	  _force.resize(newSize,Vector3r::Zero());
	  _torque.resize(newSize,Vector3r::Zero());}
  if (permForceUsed) resizePerm(newSize);
  if (moveRotUsed and _move.size()<newSize) {
	  _move.reserve(size_t(newSize*1.3));
	  _rot.reserve(size_t(newSize*1.3));
	  _move.resize(size,Vector3r::Zero());
	  _rot.resize(size,Vector3r::Zero());}
  syncedSizes=true;
  size=newSize;
}
#endif
