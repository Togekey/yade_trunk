// (c) 2018 Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>

#define YADE_MPI
#ifdef YADE_MPI

#include "Subdomain.hpp"
#include <core/Scene.hpp>
#include <core/BodyContainer.hpp>
#include <core/State.hpp>
#include <pkg/common/Sphere.hpp>
#include <core/MPIBodyContainer.hpp>
#include <lib/serialization/ObjectIO.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream_buffer.hpp>
// #include <lib/pyutil/numpy_boost.hpp>

YADE_PLUGIN((Subdomain));
CREATE_LOGGER(Subdomain);


YADE_PLUGIN((Bo1_Subdomain_Aabb)/*(Bo1_Facet_Aabb)(Bo1_Box_Aabb)*/);

void Subdomain::setMinMax()
{
	Scene* scene(Omega::instance().getScene().get());	// get scene
// 	Vector3r min, max;
	Real inf=std::numeric_limits<Real>::infinity();
	boundsMin=Vector3r(inf,inf,inf); boundsMax=Vector3r(-inf,-inf,-inf);
	if (ids.size()==0) LOG_WARN("empty subdomain!");
	if (ids.size()>0 and Body::byId(ids[0],scene)->subdomain != scene->subdomain) LOG_WARN("setMinMax executed with deprecated data (body->subdomain != scene->subdomain)");
	FOREACH(const Body::id_t& id, ids){
		const shared_ptr<Body>& b = Body::byId(id,scene);
		if(!b or !b->bound) continue;
		boundsMax=boundsMax.cwiseMax(b->bound->max);
		boundsMin=boundsMin.cwiseMin(b->bound->min);
	}
}

// inspired by Integrator::slaves_set (Integrator.hpp)
void Subdomain::intrs_set(const boost::python::list& source){
	int len=boost::python::len(source);
	intersections.clear();
	for(int i=0; i<len; i++){
		boost::python::extract<std::vector<Body::id_t> > serialGroup(source[i]);
		if (serialGroup.check()){ intersections.push_back(serialGroup()); continue; }
		cerr<<"  ... failed"<<endl;
		PyErr_SetString(PyExc_TypeError,"intersections should be provided as a list of list of ids");
		boost::python::throw_error_already_set();
	}
}

void Subdomain::mIntrs_set(const boost::python::list& source){
	int len=boost::python::len(source);
	mirrorIntersections.clear();
	for(int i=0; i<len; i++){
		boost::python::extract<std::vector<Body::id_t> > serialGroup(source[i]);
		if (serialGroup.check()){ mirrorIntersections.push_back(serialGroup()); continue; }
		cerr<<"  ... failed"<<endl;
		PyErr_SetString(PyExc_TypeError,"intersections should be provided as a list of list of ids");
		boost::python::throw_error_already_set();
	}
}

boost::python::list Subdomain::intrs_get(){
	boost::python::list ret;
	FOREACH(vector<Body::id_t >& grp, intersections){
		ret.append(boost::python::list(grp));
	}
	return ret;
}

boost::python::list Subdomain::mIntrs_get(){
	boost::python::list ret;
	FOREACH(vector<Body::id_t >& grp, mirrorIntersections){
		ret.append(boost::python::list(grp));
	}
	return ret;
}

void Bo1_Subdomain_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* b){
// 	LOG_WARN("Bo1_Subdomain_Aabb::go()")
	Subdomain* domain = static_cast<Subdomain*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb);}
	Aabb* aabb=static_cast<Aabb*>(bv.get());

	if(!scene->isPeriodic){
		aabb->min=domain->boundsMin; aabb->max=domain->boundsMax;
		return;
	} else {LOG_ERROR("to be implemented")}
	LOG_WARN("Bo1_Subdomain_Aabb::go not implemented for periodic boundaries");
}


/********************dpk********************/

std::string Subdomain::serializeMPIBodyContainer(const shared_ptr<MPIBodyContainer>& container) {
	std::string strContainer;
	boost::iostreams::back_insert_device<std::string> inserter(strContainer);
	boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
	yade::ObjectIO::save<decltype(container), boost::archive::binary_oarchive>(s, "container", container);
	s.flush();
	return strContainer;
}


shared_ptr<MPIBodyContainer> Subdomain::deSerializeMPIBodyContainer(const char* strContainer, int sizeC) {
	 shared_ptr<MPIBodyContainer> container (shared_ptr<MPIBodyContainer>  (new MPIBodyContainer()));
	 boost::iostreams::basic_array_source<char> device(strContainer, sizeC);
	 boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
	 yade::ObjectIO::load<decltype(container),boost::archive::binary_iarchive>(s, "container", container);
	 return container;
}


string Subdomain::fillContainerGetString(shared_ptr<MPIBodyContainer>& container, std::vector<Body::id_t>& ids) {
	container->insertBodyList(ids);
	std::string containerString = serializeMPIBodyContainer(container);
	return containerString;
}

void Subdomain::clearSubdomainIds(){
	ids.clear();
}

void Subdomain::setIDstoSubdomain(boost::python::list& idList ){
	unsigned int listSize = boost::python::len(idList);
	for (unsigned int i=0; i != listSize; ++i) {
	   int  b_id = boost::python::extract<int> (idList[i]);
	   ids.push_back(b_id);
	}
}

void Subdomain::getRankSize() {
	  MPI_Comm_rank(MPI_COMM_WORLD, &subdomainRank);
	  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
}

// driver function for merge operation // workers send bodies, master recieves, sets the bodies into bodycaontainer, sets interactions in interactionContainer.
void Subdomain::mergeOp() {
	sendAllBodiesToMaster();
	recvBodyContainersFromWorkers();
	if (subdomainRank==master){
	  bool setDeletedBodies = false; setBodiesToBodyContainer(recvdBodyContainers, setDeletedBodies);
		setBodyIntrs();
		bodiesSet = false; // reset flag for next merge op.
	}
}

//TODO:: add function to set interaction--> set selected interactions.

void Subdomain::setCommunicationContainers() {

    //here, we setup the serialized MPIBodyContainer (MPIBodyContainer to string). it is std::vector<std::pair<container(string), sendRank/recvRank> >

    // fill the send container based on the ids from the intersection(local) map

    //Send container
	 if (subdomainRank == master) {return ; }
	 recvRanks.clear(); sendContainer.clear();
	 unsigned int zero = 0;
	for (unsigned int i=1; i != intersections.size(); ++i){
	  if ((intersections[i].size()== zero) || (i == static_cast<unsigned int>(subdomainRank))){continue ; } // exclude self or if no intersections with others)
	  shared_ptr<MPIBodyContainer> container(shared_ptr<MPIBodyContainer> (new MPIBodyContainer()));
	  container->subdomainRank = subdomainRank;  // used to identify the origin rank at the reciever side. (maybe not needed?)
	  std::string containerString = fillContainerGetString(container, intersections[i]);
	  sendContainer.push_back(std::make_pair(containerString, i));
	}

    //Recv container, here we just need the ranks for now.
	for (unsigned int i=1; i != remoteCount.size(); ++i) {
	  if ( ( static_cast<unsigned int>(subdomainRank) == i) || (!remoteCount[i]) ){continue  ; }
	  recvRanks.push_back(i);
	}

	commContainer = true;  //flag to check if the communicationContainers are set.
}



void Subdomain::sendContainerString() {
	//send the containers.
	if (subdomainRank == master) {return ; }
	if (! commContainer ){ LOG_ERROR("communication containers are not set!"); }
	for (unsigned int i=0; i != sendContainer.size(); ++i) {
	  MPI_Request request;
	  sendString(sendContainer[i].first, sendContainer[i].second, TAG_STRING+sendContainer[i].second, request);
	  mpiReqs.push_back(request);
	}
}


void Subdomain::processContainerStrings() {
	//convert the recieved string buffers to MPIBodyContainer.
	recvdBodyContainers.clear();
	if (! containersRecvd) {LOG_ERROR ("containerStrings not recvd. Fail!"); return ; }
	for (unsigned int i=0; i != recvRanks.size(); ++i) {
	  char *cbuf = recvdCharBuff[i]; int sz = recvdStringSizes[i];
	  cbuf[sz] = '\0';
	  shared_ptr<MPIBodyContainer> cntr(deSerializeMPIBodyContainer(cbuf, sz));
	  recvdBodyContainers.push_back(cntr);
	}
	//free the pointers
	clearRecvdCharBuff(recvdCharBuff);
}

void Subdomain::sendAllBodiesToMaster() {
  // send all bodies from this subdomain to the master. Can be used for merge.
  // (note to self: this can be improved based on the bisection decomposition.)
	if(subdomainRank == master) {return;}
	shared_ptr<MPIBodyContainer> container(shared_ptr<MPIBodyContainer> (new MPIBodyContainer()));
	std::string s = fillContainerGetString(container, ids);
	sendStringBlocking(s, 0, TAG_BODY+0);
}

/********Functions exclusive to the master*************/

//TODO:: set interactions to master's interaction container from the recvd worker MPIBodyContainers.

void Subdomain::initMasterContainer(){
	if (subdomainRank != master) {return; }
	 recvdBodyContainers.resize(commSize);
	 recvRanks.resize(commSize);
	 recvdCharBuff.resize(commSize);
	 allocContainerMaster = true;
}

void Subdomain::recvBodyContainersFromWorkers() {
	if (subdomainRank != master ) {return ; }
	if (! allocContainerMaster){ initMasterContainer(); }
	for (int sourceRank=1; sourceRank != commSize; ++sourceRank){
	 int sz =  probeIncomingBlocking(sourceRank);
	 char *cbuf = new char[sz];
	  //recvStringBlocking(recvContainer[sourceRank].first, recvContainer[sourceRank].second,TAG_BODY + subdomainRank,mpiStatus[sourceRank]);
	 recvBuffBlocking(cbuf, sz , TAG_STRING+subdomainRank, sourceRank);
	 recvdCharBuff.push_back(cbuf);
	}
}

// set all body properties from the worker MPIBodyContainer

void Subdomain::setBodiesToBodyContainer(std::vector<shared_ptr<MPIBodyContainer> >& containers, bool setDeletedBodies) { // to be used when deserializing a recieved container.

	    const shared_ptr<Scene>& scene = Omega::instance().getScene();
	    shared_ptr<BodyContainer>& bodyContainer = scene->bodies;
	    for (unsigned int i=0; i != containers.size(); ++i){
	      shared_ptr<MPIBodyContainer>& mpiContainer = containers[i];
	      std::vector<shared_ptr<Body> >& bContainer = mpiContainer->bContainer;
	      for (auto bIter = bContainer.begin(); bIter != bContainer.end(); ++ bIter) {
		shared_ptr<Body> newBody = *(bIter);
		// check if the body already exists in the existing bodycontainer
		const Body::id_t& id = newBody->id;
		if ((!(*bodyContainer)[id]) &&  setDeletedBodies) {bodyContainer->insertAtId(newBody, newBody->id);} // insert to the bodycontainer
		else{ (*bodyContainer)[id] = newBody;} // copy shared pointer.}
	      }
	    }
			containers.clear(); //free the vector of MPIBodyContainers;
			bodiesSet = true;
    }




void Subdomain::setBodyIntrs() {
	if (!bodiesSet) {LOG_ERROR("MASTER PROC : Bodies are not set in Body container."); return;  }
	shared_ptr<Scene> scene= Omega::instance().getScene();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	shared_ptr<InteractionContainer>& interactionContainer = scene -> interactions;
	std::vector<shared_ptr<Body> >&  container = bodies->body; // the real bodycontainer aka std::vector<shared_ptr<Body> >
	for (auto bIter = container.begin(); bIter != container.end(); ++bIter) {
	    const shared_ptr<Body> b = *(bIter);
	    //shared_ptr<Interaction>& bIntr = b->intrs.second;
      std::map<Body::id_t, shared_ptr<Interaction> >:: iterator mIter;
      while (mIter != b->intrs.end()){
	    interactionContainer-> insert(mIter->second);
	    mIter ++;
      }
	}
}

/*********************communication functions**************/
//blocking  send and recv

void Subdomain::sendStringBlocking(std::string& s, int destRank, int tag) {

	MPI_Send(s.data(), s.size(), MPI_CHAR, destRank, tag, MPI_COMM_WORLD);
}

int Subdomain::probeIncomingBlocking(int sourceRank) {
	MPI_Status status;
	MPI_Probe(sourceRank, TAG_STRING + subdomainRank, MPI_COMM_WORLD, &status);
	int sz;
	MPI_Get_count(&status, MPI_CHAR, &sz);
	return sz;
}

void Subdomain::recvBuffBlocking(char* cbuf, int cbufSz, int tag, int sourceRank){
	MPI_Status status;
	MPI_Recv(cbuf, cbufSz, MPI_CHAR, sourceRank, tag, MPI_COMM_WORLD, &status);
}



//non-blocking calls  --> Isend, Iprobe, Irecv, Waitall();  (we will use mpi_isend + mpi_wait and then mpi_probe followed by mpi_recv)

void Subdomain::sendString(std::string& s, int destRank, int tag,  MPI_Request& request){

	//char* cbuf  = &s[0]; int len = s.size();
	int len = s.size();
	MPI_Isend(s.data(), len, MPI_CHAR,  destRank, tag , MPI_COMM_WORLD, &request);  //

}

int Subdomain::probeIncoming(int sourceRank) {
	int flag=0; MPI_Status status;
	while (! flag) {
	  MPI_Iprobe(sourceRank, TAG_STRING + subdomainRank, MPI_COMM_WORLD, &flag, &status);
	}
	int sz;
	MPI_Get_count(&status, MPI_CHAR, &sz);
	return sz; mpiStatus.push_back(status);

}

void Subdomain::recvBuff(char* cbuf, int cbufsZ, int sourceRank, MPI_Request& request) {

	MPI_Irecv(cbuf, cbufsZ, MPI_CHAR, sourceRank, TAG_STRING + subdomainRank, MPI_COMM_WORLD, &request);
}

void Subdomain::processReqs(std::vector<MPI_Request>& mpiReqs, int index) {

	for (unsigned int i = static_cast<unsigned int>(index); i != mpiReqs.size(); ++i){
	  MPI_Status status;
	  MPI_Wait(&mpiReqs[i], &status);
    }

    resetReqs(mpiReqs);
}

void Subdomain::resetReqs(std::vector<MPI_Request>& mpiReqs) {
	mpiReqs.clear();
}

void Subdomain::processReqsAll(std::vector<MPI_Request>& mpiReqs, std::vector<MPI_Status>& mpiStats) {

	for (unsigned int i = 0 ; i != mpiReqs.size(); ++i){
	  //MPI_Status status;
	  MPI_Waitall(1, &mpiReqs[i], &mpiStats[i]);
        }
        mpiStats.clear();
        resetReqs(mpiReqs);
}


void Subdomain::clearRecvdCharBuff(std::vector<char*>& rcharBuff) {
	for (std::vector<char*>::iterator cIter = rcharBuff.begin(); cIter != rcharBuff.end(); ++cIter){
	    delete (*cIter);
	  }
	rcharBuff.clear();
}


#endif
