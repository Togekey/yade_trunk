#include "yadeWrapper_02.hpp"
#include <boost/algorithm/string.hpp>

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#endif

CREATE_CPP_LOCAL_LOGGER("yadeWrapper_02.cpp");

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.
/*
Python normally iterates over object it is has __getitem__ and __len__, which BodyContainer does.
However, it will not skip removed bodies automatically, hence this iterator which does just that.
*/
pyBodyIterator::pyBodyIterator(const shared_ptr<BodyContainer>& bc)
{
	I    = bc->begin();
	Iend = bc->end();
}
pyBodyIterator   pyBodyIterator::pyIter() { return *this; }
shared_ptr<Body> pyBodyIterator::pyNext()
{
	BodyContainer::iterator ret;
	while (I != Iend) {
		ret = I;
		++I;
		if (*ret)
			return *ret;
	}
	PyErr_SetNone(PyExc_StopIteration);
	py::throw_error_already_set(); /* never reached, but makes the compiler happier */
	throw;
}

void pyBodyContainer::checkClump(shared_ptr<Body> b)
{
	if (!(b->isClump())) {
		PyErr_SetString(PyExc_TypeError, ("Error: Body" + boost::lexical_cast<string>(b->getId()) + " is not a clump.").c_str());
		py::throw_error_already_set();
	}
}

pyBodyIterator pyBodyContainer::pyIter() { return pyBodyIterator(proxee); }
pyBodyContainer::pyBodyContainer(const shared_ptr<BodyContainer>& _proxee)
        : proxee(_proxee)
{
}
// raw access to the underlying
bool pyBodyContainer::getUseRedirection(void) { return proxee->useRedirection; }
bool pyBodyContainer::getEnableRedirection(void) { return proxee->enableRedirection; }
void pyBodyContainer::setUseRedirection(bool val)
{
	if (val and not proxee->useRedirection == val)
		proxee->useRedirection = val;
	proxee->dirty = true;
	if (val)
		proxee->enableRedirection = true;
}
void pyBodyContainer::setEnableRedirection(bool val)
{
	proxee->enableRedirection = val;
	if (not val)
		proxee->useRedirection = false;
}

shared_ptr<Body> pyBodyContainer::pyGetitem(Body::id_t _id)
{
	int id = (_id >= 0 ? _id : proxee->size() + _id);
	if (id < 0 || (size_t)id >= proxee->size()) {
		PyErr_SetString(PyExc_IndexError, "Body id out of range.");
		py::throw_error_already_set(); /* make compiler happy; never reached */
		return shared_ptr<Body>();
	} else
		return (*proxee)[id];
}
Body::id_t pyBodyContainer::append(shared_ptr<Body> b)
{
	// shoud be >=0, but Body is by default created with id 0... :-|
	if (b->getId() >= 0) {
		PyErr_SetString(
		        PyExc_IndexError,
		        ("Body already has id " + boost::lexical_cast<string>(b->getId()) + " set; appending such body (for the second time) is not allowed.")
		                .c_str());
		py::throw_error_already_set();
	}
	return proxee->insert(b);
}
Body::id_t         pyBodyContainer::insertAtId(shared_ptr<Body> b, Body::id_t pos) { return proxee->insertAtId(b, pos); }
vector<Body::id_t> pyBodyContainer::appendList(vector<shared_ptr<Body>> bb)
{
	const std::lock_guard<std::mutex> lock(Omega::instance().renderMutex);
	vector<Body::id_t>                ret;
	FOREACH(shared_ptr<Body> & b, bb) { ret.push_back(append(b)); }
	return ret;
}
Body::id_t pyBodyContainer::clump(vector<Body::id_t> ids, unsigned int discretization)
{
	// create and add clump itself
	Scene*            scene(Omega::instance().getScene().get());
	shared_ptr<Body>  clumpBody = shared_ptr<Body>(new Body());
	shared_ptr<Clump> clump     = shared_ptr<Clump>(new Clump());
	clumpBody->shape            = clump;
	clumpBody->setBounded(false);
	proxee->insert(clumpBody);
	// add clump members to the clump
	FOREACH(Body::id_t id, ids)
	{
		if (Body::byId(id, scene)->isClumpMember()) {                                                 //Check, whether the body is clumpMember
			Clump::del(Body::byId(Body::byId(id, scene)->clumpId, scene), Body::byId(id, scene)); //If so, remove it from there
		}
	};

	FOREACH(Body::id_t id, ids) Clump::add(clumpBody, Body::byId(id, scene));
	Clump::updateProperties(clumpBody, discretization);
	return clumpBody->getId();
}
py::tuple pyBodyContainer::appendClump(vector<shared_ptr<Body>> bb, unsigned int discretization)
{
	// append constituent particles
	vector<Body::id_t> ids(appendList(bb));
	// clump them together (the clump fcn) and return
	return py::make_tuple(clump(ids, discretization), ids);
}
void pyBodyContainer::updateClumpProperties(py::list excludeList, unsigned int discretization)
{
	//convert excludeList to a c++ list
	vector<Body::id_t> excludeListC;
	for (int ii = 0; ii < py::len(excludeList); ii++)
		excludeListC.push_back(py::extract<Body::id_t>(excludeList[ii])());
	for (const auto& b : *proxee) {
		if (!(std::find(excludeListC.begin(), excludeListC.end(), b->getId()) != excludeListC.end())) {
			if (b->isClump())
				Clump::updateProperties(b, discretization);
		}
	}
}
void pyBodyContainer::deleteClumpMember(shared_ptr<Body> clumpBody, shared_ptr<Body> memberBody)
{ //FIXME
	const shared_ptr<Clump> clump(YADE_PTR_CAST<Clump>(clumpBody->shape));
	if (clump->members.size() == 1) {
		Clump::del(clumpBody, memberBody); //phD was not commented out
		for (unsigned i = 0; i < clump->ids.size(); i++) {
			if (clump->ids[i] == memberBody->getId()) {
				clump->ids.erase(clump->ids.begin() + i);
			}
		}
		proxee->erase(memberBody->getId(), false);
		proxee->erase(clumpBody->getId(), false);

	} else {
		Clump::del(clumpBody, memberBody); //pHD was not commented out
		for (unsigned i = 0; i < clump->ids.size(); i++) {
			if (clump->ids[i] == memberBody->getId()) {
				clump->ids.erase(clump->ids.begin() + i);
			}
		}
		Clump::updatePropertiesNonSpherical(clumpBody, /*intersecting*/ false);
		proxee->erase(memberBody->getId(), false);
	}
}
void pyBodyContainer::deleteClumpBody(shared_ptr<Body> clumpBody)
{ //FIXME
	const shared_ptr<Clump> clump(YADE_PTR_CAST<Clump>(clumpBody->shape));

	//if (clump->members.size()==0 ){
	//	proxee->erase(clumpBody->getId());
	//}else{
	Scene* scene(Omega::instance().getScene().get());
	int    totalNumber = clump->ids.size();
	int    count       = 0;
	while (count < totalNumber) {
		//for (int i=0; i<clump->ids.size(); i++){
		shared_ptr<Body> memberBody(YADE_PTR_CAST<Body>(Body::byId(clump->ids[/*i */ 0], scene)));
		deleteClumpMember(clumpBody, memberBody);
		//clump->ids.erase(clump->ids.begin()+i);
		//proxee->erase(memberBody->getId());
		count++;
	}

	proxee->erase(clumpBody->getId(), true);
	//}
}
void pyBodyContainer::addToClump(vector<Body::id_t> bids, Body::id_t cid, unsigned int discretization)
{
	Scene*           scene(Omega::instance().getScene().get()); // get scene
	shared_ptr<Body> clp = Body::byId(cid, scene);              // get clump pointer
	checkClump(clp);
	vector<Body::id_t> eraseList;
	FOREACH(Body::id_t bid, bids)
	{
		shared_ptr<Body> bp = Body::byId(bid, scene); // get body pointer
		if (bp->isClump()) {
			if (bp == clp) {
				PyErr_Warn(
				        PyExc_UserWarning,
				        ("Warning: Body " + boost::lexical_cast<string>(bid) + " and clump " + boost::lexical_cast<string>(cid)
				         + " are the same bodies. Body was not added.")
				                .c_str());
				return;
			}
			Clump::add(clp, bp); //add clump bid to clump cid
			eraseList.push_back(bid);
		} else if (bp->isClumpMember()) {
			Body::id_t       bpClumpId      = bp->clumpId;
			shared_ptr<Body> bpClumpPointer = Body::byId(bpClumpId, scene);
			if (bpClumpPointer == clp) {
				PyErr_Warn(
				        PyExc_UserWarning,
				        ("Warning: Body " + boost::lexical_cast<string>(bid) + " is already a clump member of clump "
				         + boost::lexical_cast<string>(cid) + ". Body was not added.")
				                .c_str());
				return;
			}
			Clump::add(clp, bpClumpPointer); //add clump bpClumpId to clump cid
			eraseList.push_back(bpClumpId);
		} else
			Clump::add(clp, bp); // bp must be a standalone!
	}
	Clump::updateProperties(clp, discretization);
	FOREACH(Body::id_t bid, eraseList) proxee->erase(bid, false); //erase old clumps
}
void pyBodyContainer::releaseFromClump(Body::id_t bid, Body::id_t cid, unsigned int discretization)
{
	Scene*           scene(Omega::instance().getScene().get()); // get scene
	shared_ptr<Body> bp  = Body::byId(bid, scene);              // get body pointer
	shared_ptr<Body> clp = Body::byId(cid, scene);              // get clump pointer
	checkClump(clp);
	if (bp->isClumpMember()) {
		Body::id_t bpClumpId = bp->clumpId;
		if (cid == bpClumpId) {
			const shared_ptr<Clump>&    clump   = YADE_PTR_CAST<Clump>(clp->shape);
			std::map<Body::id_t, Se3r>& members = clump->members;
			if (members.size() == 2) {
				PyErr_Warn(
				        PyExc_UserWarning,
				        ("Warning: Body " + boost::lexical_cast<string>(bid) + " not released from clump " + boost::lexical_cast<string>(cid)
				         + ", because number of clump members would get < 2!")
				                .c_str());
				return;
			}
			Clump::del(clp, bp); //release bid from cid
			Clump::updateProperties(clp, discretization);
		} else {
			PyErr_Warn(
			        PyExc_UserWarning,
			        ("Warning: Body " + boost::lexical_cast<string>(bid) + " must be a clump member of clump " + boost::lexical_cast<string>(cid)
			         + ". Body was not released.")
			                .c_str());
			return;
		}
	} else {
		PyErr_Warn(PyExc_UserWarning, ("Warning: Body " + boost::lexical_cast<string>(bid) + " is not a clump member. Body was not released.").c_str());
		return;
	}
}
py::list pyBodyContainer::replaceByClumps(py::list ctList, vector<Real> amounts, unsigned int discretization)
{
	py::list ret;
	Real     checkSum = 0.0;
	FOREACH(Real amount, amounts)
	{
		if (amount < 0.0) {
			PyErr_SetString(PyExc_ValueError, ("Error: One or more of given amounts are negative!"));
			py::throw_error_already_set();
		} else
			checkSum += amount;
	}
	if (checkSum > 1.0) {
		PyErr_SetString(
		        PyExc_ValueError, ("Error: Sum of amounts " + boost::lexical_cast<string>(checkSum) + " should not be bigger than 1.0!").c_str());
		py::throw_error_already_set();
	}
	if (static_cast<size_t>(py::len(ctList)) != static_cast<size_t>(amounts.size())) { //avoid unsigned comparison warning
		PyErr_SetString(
		        PyExc_ValueError,
		        ("Error: Length of amounts list (" + boost::lexical_cast<string>(amounts.size()) + ") differs from length of template list ("
		         + boost::lexical_cast<string>(py::len(ctList)) + ").")
		                .c_str());
		py::throw_error_already_set();
	}
	//set a random generator (code copied from pkg/dem/SpherePack.cpp):
	static boost::minstd_rand randGen((int)TimingInfo::getNow(/* get the number even if timing is disabled globally */ true));
	typedef boost::variate_generator<boost::minstd_rand&, boost::uniform_real<>> UniRandGen;
	static UniRandGen                                                            rndUnit(randGen, boost::uniform_real<>(-1, 1));

	//get number of spherical particles and a list of all spheres:
	vector<shared_ptr<Body>> sphereList;
	shared_ptr<Sphere>       sph(new Sphere);
	int                      Sph_Index = sph->getClassIndexStatic();
	for (const auto& b : *proxee)
		if ((b->shape->getClassIndex() == Sph_Index) && (b->isStandalone()))
			sphereList.push_back(b);
	int num = sphereList.size();

	//loop over templates:
	int numSphereList = num, numTemplates = amounts.size();
	for (int ii = 0; ii < numTemplates; ii++) {
		//ctList: [<ct1>,<ct2>, ...] = [<int,[double,double, ... ],[Vector3r,Vector3r, ...]>,<int,[double,double, ... ],[Vector3r,Vector3r, ...]>, ...]
		//ct: <len(relRadList),relRadList,relPosList> = <int,[double,double, ... ],[Vector3r,Vector3r, ...]> (python objects)
		//relRadList: [relRad1,relRad2, ...] (list of doubles)
		//relPosList: [relPos1,relPos2, ...] (list of vectors)

		//extract attributes from python objects:
		py::object ctTmp         = ctList[ii];
		int        numCM         = py::extract<int>(ctTmp.attr("numCM"))(); // number of clump members
		py::list   relRadListTmp = py::extract<py::list>(ctTmp.attr("relRadii"))();
		py::list   relPosListTmp = py::extract<py::list>(ctTmp.attr("relPositions"))();

		//get relative radii and positions; calculate volumes; get balance point: get axis aligned bounding box; get minimum radius;
		vector<Real>     relRadTmp(numCM), relVolTmp(numCM);
		vector<Vector3r> relPosTmp(numCM);
		Vector3r         relPosTmpMean = Vector3r::Zero();
		Real             rMin          = 1. / 0.;
		AlignedBox3r     aabb;
		for (int jj = 0; jj < numCM; jj++) {
			relRadTmp[jj] = py::extract<Real>(relRadListTmp[jj])();
			relVolTmp[jj] = (4. / 3.) * Mathr::PI * pow(relRadTmp[jj], 3.);
			relPosTmp[jj] = py::extract<Vector3r>(relPosListTmp[jj])();
			relPosTmpMean += relPosTmp[jj];
			aabb.extend(relPosTmp[jj] + Vector3r::Constant(relRadTmp[jj]));
			aabb.extend(relPosTmp[jj] - Vector3r::Constant(relRadTmp[jj]));
			rMin = min(rMin, relRadTmp[jj]);
		}
		relPosTmpMean /= numCM; //balance point

		//get volume of the clump template using regular cubic cell array inside axis aligned bounding box of the clump:
		//(some parts are duplicated from intergration algorithm in Clump::updateProperties)
		Real     dx = rMin / 5.;     //edge length of cell
		Real     dv = pow(dx, 3);    //volume of a single cell
		Vector3r x;                  //position vector (center) of cell
		Real     relVolSumTmp = 0.0; //volume of clump template
		for (x.x() = aabb.min().x() + dx / 2.; x.x() < aabb.max().x(); x.x() += dx) {
			for (x.y() = aabb.min().y() + dx / 2.; x.y() < aabb.max().y(); x.y() += dx) {
				for (x.z() = aabb.min().z() + dx / 2.; x.z() < aabb.max().z(); x.z() += dx) {
					for (int jj = 0; jj < numCM; jj++) {
						if ((x - relPosTmp[jj]).squaredNorm() < pow(relRadTmp[jj], 2)) {
							relVolSumTmp += dv;
							break;
						}
					}
				}
			}
		}

		//get pointer lists of spheres, that should be replaced:
		int                      numReplaceTmp = int(math::round(num * amounts[ii]));
		vector<shared_ptr<Body>> bpListTmp(numReplaceTmp);
		int                      a = 0, c = 0; //counters
		vector<int>              posTmp;
		FOREACH(const shared_ptr<Body>& b, sphereList)
		{
			if (c == a * numSphereList / numReplaceTmp) {
				bpListTmp[a] = b;
				a++;
				posTmp.push_back(c); //remember position in sphereList
			}
			c++;
		}
		for (int jj = 0; jj < a; jj++) {
			sphereList.erase(sphereList.begin() + posTmp[jj] - jj); //remove bodies from sphereList, that were already found
			numSphereList--;
		}

		//adapt position- and radii-informations and replace spheres from bpListTmp by clumps:
		FOREACH(const shared_ptr<Body>& b, bpListTmp)
		{
			//get sphere, that should be replaced:
			const Sphere*        sphere = YADE_CAST<Sphere*>(b->shape.get());
			shared_ptr<Material> matTmp = b->material;

			//get a random rotation quaternion:
			Quaternionr randAxisTmp = (Quaternionr)AngleAxisr(2 * Mathr::PI * rndUnit(), Vector3r(rndUnit(), rndUnit(), rndUnit()));
			randAxisTmp.normalize();

			//convert geometries in global coordinates (scaling):
			Real               scalingFactorVolume = ((4. / 3.) * Mathr::PI * pow(sphere->radius, 3.)) / relVolSumTmp;
			Real               scalingFactor1D     = pow(scalingFactorVolume, 1. / 3.); //=((vol. sphere)/(relative clump volume))^(1/3)
			vector<Vector3r>   newPosTmp(numCM);
			vector<Real>       newRadTmp(numCM);
			vector<Body::id_t> idsTmp(numCM);
			for (int jj = 0; jj < numCM; jj++) {
				newPosTmp[jj] = relPosTmp[jj] - relPosTmpMean;   //shift position, to get balance point at (0,0,0)
				newPosTmp[jj] = randAxisTmp * newPosTmp[jj];     //rotate around balance point
				newRadTmp[jj] = relRadTmp[jj] * scalingFactor1D; //scale radii
				newPosTmp[jj] = newPosTmp[jj] * scalingFactor1D; //scale position
				newPosTmp[jj] += b->state->pos;                  //translate new position to spheres center

				//create spheres:
				shared_ptr<Body> newSphere    = shared_ptr<Body>(new Body());
				newSphere->state->blockedDOFs = State::DOF_NONE;
				newSphere->state->mass        = scalingFactorVolume * relVolTmp[jj] * matTmp->density; //vol. corrected mass for clump members
				Real inertiaTmp               = 2.0 / 5.0 * newSphere->state->mass * newRadTmp[jj] * newRadTmp[jj];
				newSphere->state->inertia     = Vector3r(inertiaTmp, inertiaTmp, inertiaTmp);
				newSphere->state->pos         = newPosTmp[jj];
				newSphere->material           = matTmp;

				shared_ptr<Sphere> sphereTmp = shared_ptr<Sphere>(new Sphere());
				sphereTmp->radius            = newRadTmp[jj];
				sphereTmp->color             = Vector3r(Mathr::UnitRandom(), Mathr::UnitRandom(), Mathr::UnitRandom());
				sphereTmp->color.normalize();
				newSphere->shape = sphereTmp;

				shared_ptr<Aabb> aabbTmp = shared_ptr<Aabb>(new Aabb());
				aabbTmp->color           = Vector3r(0, 1, 0);
				newSphere->bound         = aabbTmp;
				proxee->insert(newSphere);
				LOG_DEBUG("New body (sphere) " << newSphere->id << " added.");
				idsTmp[jj] = newSphere->id;
			}
			Body::id_t newClumpId = clump(idsTmp, discretization);
			ret.append(py::make_tuple(newClumpId, idsTmp));
			erase(b->id, false);
		}
	}
	return ret;
}
Real pyBodyContainer::getRoundness(py::list excludeList)
{
	Scene*             scene(Omega::instance().getScene().get()); // get scene
	shared_ptr<Sphere> sph(new Sphere);
	int                Sph_Index = sph->getClassIndexStatic(); // get sphere index for checking if bodies are spheres
	//convert excludeList to a c++ list
	vector<Body::id_t> excludeListC;
	for (int ii = 0; ii < py::len(excludeList); ii++)
		excludeListC.push_back(py::extract<Body::id_t>(excludeList[ii])());
	Real RC_sum = 0.0; //sum of local roundnesses
	Real R1, R2, vol, dens;
	int  c = 0; //counter
	for (const auto& b : *proxee) {
		if (!(std::find(excludeListC.begin(), excludeListC.end(), b->getId()) != excludeListC.end())) {
			if ((b->shape->getClassIndex() == Sph_Index) && (b->isStandalone())) {
				RC_sum += 1.0;
				c += 1;
			}
			if (b->isClump()) {
				R2                                  = 0.0;
				dens                                = 0.0;
				vol                                 = 0.0;
				const shared_ptr<Clump>&    clump   = YADE_PTR_CAST<Clump>(b->shape);
				std::map<Body::id_t, Se3r>& members = clump->members;
				for (MemberMap::value_type& mm : members) {
					const Body::id_t&       memberId = mm.first;
					const shared_ptr<Body>& member   = Body::byId(memberId, scene);
					assert(member->isClumpMember());
					if (member->shape->getClassIndex() == Sph_Index) { //clump member should be a sphere
						const Sphere* sphere = YADE_CAST<Sphere*>(member->shape.get());
						R2
						        = max((member->state->pos - b->state->pos).norm() + sphere->radius,
						              R2); //get minimum radius of a sphere, that imbeds clump
						dens = member->material->density;
					}
				}
				if (dens > 0.)
					vol = b->state->mass / dens;
				R1 = pow((3. * vol) / (4. * Mathr::PI), 1. / 3.); //get theoretical radius of a sphere, with same volume as clump
				if (R2 < R1) {
					PyErr_Warn(PyExc_UserWarning, ("Something went wrong in getRoundness method (R2 < R1 detected)."));
					return 0;
				}
				RC_sum += R1 / R2;
				c += 1;
			}
		}
	}
	if (c == 0)
		c = 1;     //in case no spheres and no clumps are present in the scene: RC = 0
	return RC_sum / c; //return roundness coefficient RC
}
vector<Body::id_t> pyBodyContainer::replace(vector<shared_ptr<Body>> bb)
{
	proxee->clear();
	return appendList(bb);
}
long pyBodyContainer::length() { return proxee->size(); }
void pyBodyContainer::clear() { proxee->clear(); }
bool pyBodyContainer::erase(Body::id_t id, bool eraseClumpMembers) { return proxee->erase(id, eraseClumpMembers); }
#ifdef YADE_MPI
vector<Body::id_t> pyBodyContainer::subdomainBodies() { return proxee->subdomainBodies; }
#endif


pyTags::pyTags(const shared_ptr<Scene> _mb)
        : mb(_mb)
{
}
bool pyTags::hasKey(const string& key)
{
	FOREACH(string val, mb->tags)
	{
		if (boost::algorithm::starts_with(val, key + "=")) {
			return true;
		}
	}
	return false;
}
string pyTags::getItem(const string& key)
{
	FOREACH(string & val, mb->tags)
	{
		if (boost::algorithm::starts_with(val, key + "=")) {
			string val1(val);
			boost::algorithm::erase_head(val1, key.size() + 1);
			return val1;
		}
	}
	PyErr_SetString(PyExc_KeyError, ("Invalid key: " + key + ".").c_str());
	py::throw_error_already_set(); /* make compiler happy; never reached */
	return string();
}
void pyTags::setItem(const string& key, const string& item)
{
	if (key.find("=") != string::npos) {
		PyErr_SetString(PyExc_KeyError, "Key must not contain the '=' character (implementation limitation; sorry).");
		py::throw_error_already_set();
	}
	FOREACH(string & val, mb->tags)
	{
		if (boost::algorithm::starts_with(val, key + "=")) {
			val = key + "=" + item;
			return;
		}
	}
	mb->tags.push_back(key + "=" + item);
}
py::list pyTags::keys()
{
	py::list ret;
	FOREACH(string val, mb->tags)
	{
		size_t i = val.find("=");
		if (i == string::npos)
			throw runtime_error("Tags must be in the key=value format (internal error?)");
		boost::algorithm::erase_tail(val, val.size() - i);
		ret.append(val);
	}
	return ret;
}


pyInteractionIterator::pyInteractionIterator(const shared_ptr<InteractionContainer>& ic)
{
	I    = ic->begin();
	Iend = ic->end();
}
pyInteractionIterator   pyInteractionIterator::pyIter() { return *this; }
shared_ptr<Interaction> pyInteractionIterator::pyNext()
{
	InteractionContainer::iterator ret;
	while (I != Iend) {
		ret = I;
		++I;
		if ((*ret)->isReal())
			return *ret;
	}
	PyErr_SetNone(PyExc_StopIteration);
	py::throw_error_already_set();
	throw;  // to avoid compiler warning; never reached
	        // InteractionContainer::iterator ret=I; ++I; return *ret;
}

pyInteractionContainer::pyInteractionContainer(const shared_ptr<InteractionContainer>& _proxee)
        : proxee(_proxee)
        , scene(Omega::instance().getScene())
{
}
pyInteractionIterator   pyInteractionContainer::pyIter() { return pyInteractionIterator(proxee); }
bool                    pyInteractionContainer::has(Body::id_t id1, Body::id_t id2) { return proxee->found(id1, id2); }
shared_ptr<Interaction> pyInteractionContainer::pyGetitem(vector<Body::id_t> id12)
{
	//if(!PySequence_Check(id12.ptr())) throw invalid_argument("Key must be a tuple");
	//if(py::len(id12)!=2) throw invalid_argument("Key must be a 2-tuple: id1,id2.");
	if (id12.size() == 2) {
		//if(max(id12[0],id12[1])>
		shared_ptr<Interaction> i = proxee->find(id12[0], id12[1]);
		if (i)
			return i;
		else {
			PyErr_SetString(PyExc_IndexError, "No such interaction");
			py::throw_error_already_set(); /* make compiler happy; never reached */
			return shared_ptr<Interaction>();
		}
	} else if (id12.size() == 1) {
		return (*proxee)[id12[0]];
	} else
		throw invalid_argument("2 integers (id1,id2) or 1 integer (nth) required.");
}
/* return nth _real_ iteration from the container (0-based index); this is to facilitate picking random interaction */
shared_ptr<Interaction> pyInteractionContainer::pyNth(long n)
{
	long i = 0;
	FOREACH(shared_ptr<Interaction> I, *proxee)
	{
		if (!I->isReal())
			continue;
		if (i++ == n)
			return I;
	}
	PyErr_SetString(
	        PyExc_IndexError,
	        (string("Interaction number out of range (") + boost::lexical_cast<string>(n) + ">=" + boost::lexical_cast<string>(i) + ").").c_str());
	py::throw_error_already_set(); /* make compiler happy; never reached */
	return shared_ptr<Interaction>();
}
long     pyInteractionContainer::len() { return proxee->size(); }
void     pyInteractionContainer::clear() { proxee->clear(); }
py::list pyInteractionContainer::withBody(long id)
{
	py::list ret;
	FOREACH(const Body::MapId2IntrT::value_type& I, Body::byId(id, scene)->intrs)
	{
		if (I.second->isReal())
			ret.append(I.second);
	}
	return ret;
}
py::list pyInteractionContainer::withBodyAll(long id)
{
	py::list ret;
	FOREACH(const Body::MapId2IntrT::value_type& I, Body::byId(id, scene)->intrs) ret.append(I.second);
	return ret;
}
py::list pyInteractionContainer::getAll(bool onlyReal)
{
	py::list ret;
	FOREACH(const shared_ptr<Interaction>& I, *proxee)
	{
		if (onlyReal && !I->isReal())
			continue;
		ret.append(I);
	}
	return ret;
}
long pyInteractionContainer::countReal()
{
	long ret = 0;
	FOREACH(const shared_ptr<Interaction>& I, *proxee)
	{
		if (I->isReal())
			ret++;
	}
	return ret;
}
bool pyInteractionContainer::serializeSorted_get() { return proxee->serializeSorted; }
void pyInteractionContainer::serializeSorted_set(bool ss) { proxee->serializeSorted = ss; }
void pyInteractionContainer::eraseNonReal() { proxee->eraseNonReal(); }
void pyInteractionContainer::erase(Body::id_t id1, Body::id_t id2) { proxee->requestErase(id1, id2); }

pyForceContainer::pyForceContainer(shared_ptr<Scene> _scene)
        : scene(_scene)
{
}
void pyForceContainer::checkId(long id)
{
	if (id < 0 || (size_t)id >= scene->bodies->size()) {
		PyErr_SetString(PyExc_IndexError, "Body id out of range.");
		py::throw_error_already_set(); /* never reached */
		throw;
	}
}
Vector3r pyForceContainer::force_get(long id, bool sync)
{
	checkId(id);
	if (!sync and !scene->forces.synced)
		return scene->forces.getForceSingle(id);
	scene->forces.sync();
	return scene->forces.getForce(id);
}
Vector3r pyForceContainer::torque_get(long id, bool sync)
{
	checkId(id);
	if (!sync and !scene->forces.synced)
		return scene->forces.getTorqueSingle(id);
	scene->forces.sync();
	return scene->forces.getTorque(id);
}
void pyForceContainer::force_add(long id, const Vector3r& f, bool permanent)
{
	checkId(id);
	if (!permanent)
		scene->forces.addForce(id, f);
	else {
		LOG_WARN("O.forces.addF(...,permanent=True) is deprecated, use O.forces.setPermF(...) instead");
		scene->forces.setPermForce(id, f);
	}
}
void pyForceContainer::torque_add(long id, const Vector3r& t, bool permanent)
{
	checkId(id);
	if (!permanent)
		scene->forces.addTorque(id, t);
	else {
		LOG_WARN("O.forces.addT(...,permanent=True) is deprecated, use O.forces.setPermT(...) instead");
		scene->forces.setPermTorque(id, t);
	}
}
Vector3r pyForceContainer::permForce_get(long id)
{
	checkId(id);
	return scene->forces.getPermForce(id);
}
Vector3r pyForceContainer::permTorque_get(long id)
{
	checkId(id);
	return scene->forces.getPermTorque(id);
}
void pyForceContainer::permForce_set(long id, const Vector3r& f)
{
	checkId(id);
	scene->forces.setPermForce(id, f);
}
void pyForceContainer::permTorque_set(long id, const Vector3r& t)
{
	checkId(id);
	scene->forces.setPermTorque(id, t);
}
void pyForceContainer::reset(bool resetAll) { scene->forces.reset(scene->iter, resetAll); }
long pyForceContainer::syncCount_get() { return scene->forces.syncCount; }
void pyForceContainer::syncCount_set(long count) { scene->forces.syncCount = count; }
bool pyForceContainer::getPermForceUsed() { return scene->forces.getPermForceUsed(); }

pyMaterialContainer::pyMaterialContainer(shared_ptr<Scene> _scene)
        : scene(_scene)
{
}
shared_ptr<Material> pyMaterialContainer::getitem_id(int _id)
{
	int id = (_id >= 0 ? _id : scene->materials.size() + _id);
	if (id < 0 || (size_t)id >= scene->materials.size()) {
		PyErr_SetString(PyExc_IndexError, "Material id out of range.");
		py::throw_error_already_set(); /* never reached */
		throw;
	}
	return Material::byId(id, scene);
}
shared_ptr<Material> pyMaterialContainer::getitem_label(string label)
{
	// translate runtime_error to KeyError (instead of RuntimeError) if the material doesn't exist
	try {
		return Material::byLabel(label, scene);
	} catch (std::runtime_error& e) {
		PyErr_SetString(PyExc_KeyError, e.what());
		py::throw_error_already_set(); /* never reached; avoids warning */
		throw;
	}
}
int pyMaterialContainer::append(shared_ptr<Material> m)
{
	scene->materials.push_back(m);
	m->id = scene->materials.size() - 1;
	return m->id;
}
vector<int> pyMaterialContainer::appendList(vector<shared_ptr<Material>> mm)
{
	vector<int> ret;
	FOREACH(shared_ptr<Material> & m, mm) ret.push_back(append(m));
	return ret;
}
int pyMaterialContainer::len() { return (int)scene->materials.size(); }
int pyMaterialContainer::index(const std::string& label) { return Material::byLabelIndex(label, scene.get()); }

} // namespace yade

