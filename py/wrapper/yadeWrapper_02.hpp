#include <core/BodyContainer.hpp>
#include <core/Clump.hpp>
#include <core/InteractionContainer.hpp>
#include <pkg/common/Sphere.hpp>

#pragma once
namespace yade { // Cannot have #include directive inside.

class pyBodyIterator {
	BodyContainer::iterator I, Iend;

public:
	pyBodyIterator(const shared_ptr<BodyContainer>& bc);
	pyBodyIterator   pyIter();
	shared_ptr<Body> pyNext();
};

class pyBodyContainer {
private:
	void                               checkClump(shared_ptr<Body> b);
	typedef std::map<Body::id_t, Se3r> MemberMap;

public:
	const shared_ptr<BodyContainer> proxee;
	pyBodyIterator                  pyIter();
	pyBodyContainer(const shared_ptr<BodyContainer>& _proxee);
	// raw access to the underlying
	bool               getUseRedirection(void);
	bool               getEnableRedirection(void);
	void               setUseRedirection(bool val);
	void               setEnableRedirection(bool val);
	shared_ptr<Body>   pyGetitem(Body::id_t _id);
	Body::id_t         append(shared_ptr<Body> b);
	Body::id_t         insertAtId(shared_ptr<Body> b, Body::id_t pos);
	vector<Body::id_t> appendList(vector<shared_ptr<Body>> bb);
	Body::id_t         clump(vector<Body::id_t> ids, unsigned int discretization);
	py::tuple          appendClump(vector<shared_ptr<Body>> bb, unsigned int discretization);
	void               updateClumpProperties(py::list excludeList, unsigned int discretization);
	void               deleteClumpMember(shared_ptr<Body> clumpBody, shared_ptr<Body> memberBody);
	void               deleteClumpBody(shared_ptr<Body> clumpBody);
	void               addToClump(vector<Body::id_t> bids, Body::id_t cid, unsigned int discretization);
	void               releaseFromClump(Body::id_t bid, Body::id_t cid, unsigned int discretization);
	py::list           replaceByClumps(py::list ctList, vector<Real> amounts, unsigned int discretization);
	Real               getRoundness(py::list excludeList);
	vector<Body::id_t> replace(vector<shared_ptr<Body>> bb);
	long               length();
	void               clear();
	bool               erase(Body::id_t id, bool eraseClumpMembers);
#ifdef YADE_MPI
	vector<Body::id_t> subdomainBodies();
#endif
};

class pyTags {
public:
	pyTags(const shared_ptr<Scene> _mb);
	const shared_ptr<Scene> mb;
	bool                    hasKey(const string& key);
	string                  getItem(const string& key);
	void                    setItem(const string& key, const string& item);
	py::list                keys();
};


class pyInteractionIterator {
	InteractionContainer::iterator I, Iend;

public:
	pyInteractionIterator(const shared_ptr<InteractionContainer>& ic);
	pyInteractionIterator   pyIter();
	shared_ptr<Interaction> pyNext();
};

class pyInteractionContainer {
public:
	const shared_ptr<InteractionContainer> proxee;
	const shared_ptr<Scene>                scene;
	pyInteractionContainer(const shared_ptr<InteractionContainer>& _proxee);
	pyInteractionIterator   pyIter();
	bool                    has(Body::id_t id1, Body::id_t id2);
	shared_ptr<Interaction> pyGetitem(vector<Body::id_t> id12);
	/* return nth _real_ iteration from the container (0-based index); this is to facilitate picking random interaction */
	shared_ptr<Interaction> pyNth(long n);
	long                    len();
	void                    clear();
	py::list                withBody(long id);
	py::list                withBodyAll(long id);
	py::list                getAll(bool onlyReal);
	long                    countReal();
	bool                    serializeSorted_get();
	void                    serializeSorted_set(bool ss);
	void                    eraseNonReal();
	void                    erase(Body::id_t id1, Body::id_t id2);
};

class pyForceContainer {
	shared_ptr<Scene> scene;

public:
	pyForceContainer(shared_ptr<Scene> _scene);
	void     checkId(long id);
	Vector3r force_get(long id, bool sync);
	Vector3r torque_get(long id, bool sync);
	void     force_add(long id, const Vector3r& f, bool permanent);
	void     torque_add(long id, const Vector3r& t, bool permanent);
	Vector3r permForce_get(long id);
	Vector3r permTorque_get(long id);
	void     permForce_set(long id, const Vector3r& f);
	void     permTorque_set(long id, const Vector3r& t);
	void     reset(bool resetAll);
	long     syncCount_get();
	void     syncCount_set(long count);
	bool     getPermForceUsed();
};

class pyMaterialContainer {
	shared_ptr<Scene> scene;

public:
	pyMaterialContainer(shared_ptr<Scene> _scene);
	shared_ptr<Material> getitem_id(int _id);
	shared_ptr<Material> getitem_label(string label);
	int                  append(shared_ptr<Material> m);
	vector<int>          appendList(vector<shared_ptr<Material>> mm);
	int                  len();
	int                  index(const std::string& label);
};

} // namespace yade

