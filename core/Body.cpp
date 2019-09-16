
#include<core/Body.hpp>
#include<core/Scene.hpp>
#include<core/Omega.hpp>
#include<core/InteractionContainer.hpp>
#include<pkg/common/InsertionSortCollider.hpp>

CREATE_LOGGER(Body);

//! This could be -1 if id_t is re-typedef'ed as `int'
const Body::id_t Body::ID_NONE=Body::id_t(-1);
bool Body::shortListCheckedOnce=false;

const shared_ptr<Body>& Body::byId(Body::id_t _id, Scene* rb){return (*((rb?rb:Omega::instance().getScene().get())->bodies))[_id];}
const shared_ptr<Body>& Body::byId(Body::id_t _id, shared_ptr<Scene> rb){return (*(rb->bodies))[_id];}

// return list of interactions of this particle
boost::python::list Body::py_intrs(){
  boost::python::list ret;
	for(Body::MapId2IntrT::iterator it=this->intrs.begin(),end=this->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		ret.append((*it).second);
	}
	return ret;
}

// return number of interactions of this particle
unsigned int Body::coordNumber() const {
	unsigned int intrSize = 0;
	for(auto it=this->intrs.begin(),end=this->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		intrSize++;
	}
	return intrSize;
}

void Body::setBounded(bool d) {
	if(d) flags|=FLAG_BOUNDED; else flags&=~(FLAG_BOUNDED);
// 	if (id!=Body::ID_NONE and not shortListCheckedOnce) {
// 		shared_ptr<InsertionSortCollider> isc;
// 		FOREACH(shared_ptr<Engine>& e, Omega::instance().getScene()->engines){ isc=YADE_PTR_DYN_CAST<InsertionSortCollider>(e); if(isc) break; }
// 		if(isc and isc->keepListsShort) LOG_ERROR("changing body::bounded after insertion is not supported if collider::keepListsShort=True, turn keepListsShort off (suboptimal)");
// 		shortListCheckedOnce=true;
// 	}
}

bool Body::maskOk(int mask) const { return (mask==0 || ((groupMask & mask) != 0)); }
bool Body::maskCompatible(int mask) const { return (groupMask & mask) != 0; }
#ifdef YADE_MASK_ARBITRARY
bool Body::maskOk(const mask_t& mask) const { return (mask==0 || ((groupMask & mask) != 0)); }
bool Body::maskCompatible(const mask_t& mask) const { return (groupMask & mask) != 0; }
#endif



