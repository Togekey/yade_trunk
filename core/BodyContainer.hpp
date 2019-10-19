// 2004 © Olivier Galizzi <olivier.galizzi@imag.fr>
// 2004 © Janek Kozicki <cosurgi@berlios.de>
// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include <boost/tuple/tuple.hpp>
#include <lib/serialization/Serializable.hpp>

namespace yade { // Cannot have #include directive inside.

class Body;
class InteractionContainer;

#ifdef YADE_OPENMP
#define YADE_PARALLEL_FOREACH_BODY_BEGIN(b_, bodies)                                                                                       \
	bodies->updateShortLists();                                                                                                        \
	const vector<Body::id_t>& realBodies = bodies->realBodies;                                                                         \
	const bool redirect = bodies->useRedirection;                                                                                      \
	const Body::id_t _sz(redirect ? realBodies.size() : bodies->size());                                                               \
	_Pragma("omp parallel for") for (int k = 0; k < _sz; k++)                                                                          \
	{                                                                                                                                  \
		if (not redirect and not(*bodies)[k])                                                                                      \
			continue;                                                                                                          \
		b_((*bodies)[redirect ? realBodies[k] : k]);
#else
#define YADE_PARALLEL_FOREACH_BODY_BEGIN(b_, bodies)                                                                                       \
	bodies->updateShortLists();                                                                                                        \
	const vector<Body::id_t>& realBodies = bodies->realBodies;                                                                         \
	const bool redirect = bodies->useRedirection;                                                                                      \
	const Body::id_t _sz(redirect ? realBodies.size() : bodies->size());                                                               \
	for (int k = 0; k < _sz; k++) {                                                                                                    \
		if (not redirect and not(*bodies)[k])                                                                                      \
			continue;                                                                                                          \
		b_((*bodies)[redirect ? realBodies[k] : k]);
#endif
#define YADE_PARALLEL_FOREACH_BODY_END() }

/*
Container of bodies implemented as flat std::vector.
The nested iterators and the specialized FOREACH_BODY macros above will silently skip null body pointers which may exist after removal. The
null pointers can still be accessed via the [] operator.

Any alternative implementation should use the same API.
*/
class BodyContainer : public Serializable {
    private:
	using ContainerT = std::vector<shared_ptr<Body>>;

    public:
	friend class InteractionContainer; // accesses the body vector directly

	Body::id_t insert(shared_ptr<Body>); // => body.push_back()
	Body::id_t insertAtId(shared_ptr<Body> b, Body::id_t candidate); // => body[candidate]=...

	void clear();
	auto begin() const { return body.begin(); }
	auto end() const { return body.end(); }

	size_t size() const { return body.size(); }
	const auto & operator[](unsigned int id) const { return body[id]; }

	bool exists(Body::id_t id) const { return ((id >= 0) && ((size_t)id < body.size()) && ((bool)body[id])); }
	bool erase(Body::id_t id, bool eraseClumpMembers);

	void updateShortLists();

	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BodyContainer, Serializable, "Standard body container for a scene",
	    ((ContainerT, body, , , "The underlying vector<shared_ptr<Body> >"))(
	        (bool, dirty, true, (Attr::noSave | Attr::readonly | Attr::hidden),
	            "true if after insertion/removal of bodies, used only if collider::keepListsShort=true"))(
	        (bool, checkedByCollider, false, (Attr::noSave | Attr::readonly | Attr::hidden), ""))((vector<Body::id_t>, insertedBodies,
	        vector<Body::id_t>(), Attr::readonly, "The list of newly bodies inserted, to be used and purged by collider"))(
	        (vector<Body::id_t>, realBodies, vector<Body::id_t>(), Attr::readonly,
	            "Redirection vector to non-null bodies, used to optimize loops after numerous insertion/erase. In MPI runs the list is "
	            "restricted to bodies and neighbors present in current subdomain."))((bool, useRedirection, false, ,
	        "true if the scene uses up-to-date lists for boundedBodies and realBodies; turned true automatically 1/ after removal of "
	        "bodies if :yref:`enableRedirection`=True, and 2/ in MPI execution."))((bool, enableRedirection, true, ,
	        "let collider switch to optimized algorithm with body redirection when bodies are erased - true by default"))
#ifdef YADE_MPI
	        ((vector<Body::id_t>, subdomainBodies, vector<Body::id_t>(), , "The list of bounded bodies in the subdomain"))
#endif
	            ,
	    /*ctor*/,
	    .def("updateShortLists", &BodyContainer::updateShortLists,
	        "update lists realBodies and subdomainBodies. This function is called automatically by e.g. ForceContainer::reset(), it is "
	        "safe to call multiple times from many places since if the lists are up-to-date he function will just return."))

	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(BodyContainer);

} // namespace yade
