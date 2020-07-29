// 2007,2008 © Václav Šmilauer <eudoxos@arcig.cz>

#include <lib/base/AliasNamespaces.hpp>
#include <lib/base/Logging.hpp>
#include <lib/pyutil/doc_opts.hpp>
#include <lib/pyutil/gil.hpp>
#include <lib/pyutil/raw_constructor.hpp>
#include <lib/serialization/ObjectIO.hpp>
#include <core/Clump.hpp>
#include <core/EnergyTracker.hpp>
#include <core/FileGenerator.hpp>
#include <core/Functor.hpp>
#include <core/GlobalEngine.hpp>
#include <core/Omega.hpp>
#include <core/PartialEngine.hpp>
#include <core/ThreadRunner.hpp>
#include <core/Timing.hpp>
#include <pkg/common/Collider.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/common/InteractionLoop.hpp>
#include <pkg/common/KinematicEngines.hpp>
#include <pkg/common/ParallelEngine.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/LubricationWithPotential.hpp>
#include <pkg/dem/STLImporter.hpp>
#include <boost/archive/codecvt_null.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <csignal>
#include <locale>

#include "yadeWrapper_02.hpp"

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#endif

CREATE_CPP_LOCAL_LOGGER("yadeWrapper.cpp");

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

void termHandlerNormal(int /*sig*/)
{
	cerr << "Yade: normal exit." << endl;
	raise(SIGTERM);
}
void termHandlerError(int /*sig*/)
{
	cerr << "Yade: error exit." << endl;
	raise(SIGTERM);
}

class pyOmega {
private:
	// can be safely removed now, since pyOmega makes an empty scene in the constructor, if there is none
	void assertScene()
	{
		if (!OMEGA.getScene())
			throw std::runtime_error("No Scene instance?!");
	}
	Omega& OMEGA;

public:
	pyOmega()
	        : OMEGA(Omega::instance())
	{
		shared_ptr<Scene> rb = OMEGA.getScene();
		if (!rb) {
			OMEGA.init();
			rb = OMEGA.getScene();
		}
		assert(rb);
		if (!OMEGA.hasSimulationLoop()) {
			OMEGA.createSimulationLoop();
		}
	};
	/* Create variables in python's __builtin__ namespace that correspond to labeled objects. At this moment, only engines and functors can be labeled (not bodies etc). */
	void mapLabeledEntitiesToVariables()
	{
// not sure if we should map materials to variables by default...
// a call to this functions would have to be added to pyMaterialContainer::append
#if 0
			FOREACH(const shared_ptr<Material>& m, OMEGA.getScene()->materials){
				if(!m->label.empty()) { PyGILState_STATE gstate; gstate = PyGILState_Ensure(); PyRun_SimpleString(("__builtins__."+m->label+"=Omega().materials["+boost::lexical_cast<string>(m->id)+"]").c_str()); PyGILState_Release(gstate); }
			}
#endif
		FOREACH(const shared_ptr<Engine>& e, OMEGA.getScene()->engines)
		{
			if (!e->label.empty()) {
				pyRunString("__builtins__." + e->label + "=Omega().labeledEngine('" + e->label + "')");
			}
#define _DO_FUNCTORS(functors, FunctorT)                                                                                                                       \
	{                                                                                                                                                      \
		FOREACH(const shared_ptr<FunctorT>& f, functors)                                                                                               \
		{                                                                                                                                              \
			if (!f->label.empty()) {                                                                                                               \
				pyRunString("__builtins__." + f->label + "=Omega().labeledEngine('" + f->label + "')");                                        \
			}                                                                                                                                      \
		}                                                                                                                                              \
	}
#define _TRY_DISPATCHER(DispatcherT)                                                                                                                           \
	{                                                                                                                                                      \
		DispatcherT* d = dynamic_cast<DispatcherT*>(e.get());                                                                                          \
		if (d) {                                                                                                                                       \
			_DO_FUNCTORS(d->functors, DispatcherT::FunctorType);                                                                                   \
		}                                                                                                                                              \
	}
			_TRY_DISPATCHER(BoundDispatcher);
			_TRY_DISPATCHER(IGeomDispatcher);
			_TRY_DISPATCHER(IPhysDispatcher);
			_TRY_DISPATCHER(LawDispatcher);
			InteractionLoop* id = dynamic_cast<InteractionLoop*>(e.get());
			if (id) {
				_DO_FUNCTORS(id->geomDispatcher->functors, IGeomFunctor);
				_DO_FUNCTORS(id->physDispatcher->functors, IPhysFunctor);
				_DO_FUNCTORS(id->lawDispatcher->functors, LawFunctor);
			}
			Collider* coll = dynamic_cast<Collider*>(e.get());
			if (coll) {
				_DO_FUNCTORS(coll->boundDispatcher->functors, BoundFunctor);
			}
#undef _DO_FUNCTORS
#undef _TRY_DISPATCHER
			CombinedKinematicEngine* cke = dynamic_cast<CombinedKinematicEngine*>(e.get());
			if (cke) {
				FOREACH(const shared_ptr<KinematicEngine>& ke, cke->comb)
				{
					if (!ke->label.empty()) {
						pyRunString("__builtins__." + ke->label + "=Omega().labeledEngine('" + ke->label + "')");
					}
				}
			}
		}
	}
	py::object labeled_engine_get(string label)
	{
		FOREACH(const shared_ptr<Engine>& e, OMEGA.getScene()->engines)
		{
#define _DO_FUNCTORS(functors, FunctorT)                                                                                                                       \
	{                                                                                                                                                      \
		FOREACH(const shared_ptr<FunctorT>& f, functors)                                                                                               \
		{                                                                                                                                              \
			if (f->label == label)                                                                                                                 \
				return py::object(f);                                                                                                          \
		}                                                                                                                                              \
	}
#define _TRY_DISPATCHER(DispatcherT)                                                                                                                           \
	{                                                                                                                                                      \
		DispatcherT* d = dynamic_cast<DispatcherT*>(e.get());                                                                                          \
		if (d) {                                                                                                                                       \
			_DO_FUNCTORS(d->functors, DispatcherT::FunctorType);                                                                                   \
		}                                                                                                                                              \
	}
			if (e->label == label) {
				return py::object(e);
			}
			_TRY_DISPATCHER(BoundDispatcher);
			_TRY_DISPATCHER(IGeomDispatcher);
			_TRY_DISPATCHER(IPhysDispatcher);
			_TRY_DISPATCHER(LawDispatcher);
			InteractionLoop* id = dynamic_cast<InteractionLoop*>(e.get());
			if (id) {
				_DO_FUNCTORS(id->geomDispatcher->functors, IGeomFunctor);
				_DO_FUNCTORS(id->physDispatcher->functors, IPhysFunctor);
				_DO_FUNCTORS(id->lawDispatcher->functors, LawFunctor);
			}
			Collider* coll = dynamic_cast<Collider*>(e.get());
			if (coll) {
				_DO_FUNCTORS(coll->boundDispatcher->functors, BoundFunctor);
			}
#undef _DO_FUNCTORS
#undef _TRY_DISPATCHER
			CombinedKinematicEngine* cke = dynamic_cast<CombinedKinematicEngine*>(e.get());
			if (cke) {
				FOREACH(const shared_ptr<KinematicEngine>& ke, cke->comb)
				{
					if (ke->label == label)
						return py::object(ke);
				}
			}
		}
		throw std::invalid_argument(string("No engine labeled `") + label + "'");
	}

	long iter() { return OMEGA.getScene()->iter; }
	int  subStep() { return OMEGA.getScene()->subStep; }
	bool subStepping_get() { return OMEGA.getScene()->subStepping; }
	void subStepping_set(bool val) { OMEGA.getScene()->subStepping = val; }

	Real time() { return OMEGA.getScene()->time; }
	Real realTime() { return OMEGA.getRealTime(); }
	Real speed() { return OMEGA.getScene()->speed; }
	Real dt_get() { return OMEGA.getScene()->dt; }
	void dt_set(Real dt)
	{
		Scene* scene = OMEGA.getScene().get();
		// activate timestepper, if possible (throw exception if there is none)
		if (dt < 0) {
			if (!scene->timeStepperActivate(true)) /* not activated*/
				throw runtime_error("No TimeStepper found in O.engines.");
		} else {
			scene->dt = dt;
		}
	}
	bool dynDt_get() { return OMEGA.getScene()->timeStepperActive(); }
	bool dynDt_set(bool activate)
	{
		if (not OMEGA.getScene()->timeStepperActivate(activate) and activate) /* not activated despite passed value */
			throw runtime_error("No TimeStepper found in O.engines.");
		return true;
	}
	bool dynDtAvailable_get() { return OMEGA.getScene()->timeStepperPresent(); }
	long stopAtIter_get() { return OMEGA.getScene()->stopAtIter; }
	void stopAtIter_set(long s) { OMEGA.getScene()->stopAtIter = s; }
	Real stopAtTime_get() { return OMEGA.getScene()->stopAtTime; }
	void stopAtTime_set(Real s) { OMEGA.getScene()->stopAtTime = s; }


	bool timingEnabled_get() { return TimingInfo::enabled; }
	void timingEnabled_set(bool enabled) { TimingInfo::enabled = enabled; }
	// deprecated:
	unsigned long forceSyncCount_get() { return OMEGA.getScene()->forces.syncCount; }
	void          forceSyncCount_set(unsigned long count) { OMEGA.getScene()->forces.syncCount = count; }

	void run(long int numIter = -1, bool doWait = false)
	{
		Scene* scene = OMEGA.getScene().get();
		if (numIter > 0)
			scene->stopAtIter = scene->iter + numIter;
		OMEGA.run();
		// timespec t1,t2; t1.tv_sec=0; t1.tv_nsec=40000000; /* 40 ms */
		// while(!OMEGA.isRunning()) nanosleep(&t1,&t2); // wait till we start, so that calling wait() immediately afterwards doesn't return immediately
		LOG_DEBUG(
		        "RUN"
		        << ((scene->stopAtIter - scene->iter) > 0 ? string(" (" + boost::lexical_cast<string>(scene->stopAtIter - scene->iter) + " to go)")
		                                                  : string(""))
		        << "!");
		if (doWait)
			wait();
	}
	void pause()
	{
		Py_BEGIN_ALLOW_THREADS;
		OMEGA.pause();
		Py_END_ALLOW_THREADS;
		LOG_DEBUG("PAUSE!");
	}
	void step()
	{
		if (OMEGA.isRunning())
			throw runtime_error("Called O.step() while simulation is running.");
		OMEGA.getScene()->moveToNextTimeStep(); /* LOG_DEBUG("STEP!"); run(1); wait(); */
	}
	void wait()
	{
		if (OMEGA.isRunning()) {
			LOG_DEBUG("WAIT!");
		} else
			return;
		timespec t1, t2;
		t1.tv_sec  = 0;
		t1.tv_nsec = 40000000; /* 40 ms */
		Py_BEGIN_ALLOW_THREADS;
		while (OMEGA.isRunning())
			nanosleep(&t1, &t2);
		Py_END_ALLOW_THREADS;
		if (!OMEGA.simulationLoop->workerThrew)
			return;
		LOG_ERROR("Simulation error encountered.");
		OMEGA.simulationLoop->workerThrew = false;
		throw OMEGA.simulationLoop->workerException;
	}
	bool       isRunning() { return OMEGA.isRunning(); }
	py::object get_filename()
	{
		string f = OMEGA.sceneFile;
		if (f.size() > 0)
			return py::object(f);
		return py::object();
	}
	void load(std::string fileName, bool quiet = false)
	{
		Py_BEGIN_ALLOW_THREADS;
		OMEGA.stop();
		Py_END_ALLOW_THREADS;
		OMEGA.loadSimulation(fileName, quiet);
		OMEGA.createSimulationLoop();
		mapLabeledEntitiesToVariables();
	}
	void     reload(bool quiet = false) { load(OMEGA.sceneFile, quiet); }
	void     saveTmp(string mark = "", bool quiet = false) { save(":memory:" + mark, quiet); }
	void     loadTmp(string mark = "", bool quiet = false) { load(":memory:" + mark, quiet); }
	py::list lsTmp()
	{
		py::list                          ret;
		typedef pair<std::string, string> strstr;
		FOREACH(const strstr& sim, OMEGA.memSavedSimulations)
		{
			string mark = sim.first;
			boost::algorithm::replace_first(mark, ":memory:", "");
			ret.append(mark);
		}
		return ret;
	}
	void tmpToFile(string mark, string filename)
	{
		if (OMEGA.memSavedSimulations.count(":memory:" + mark) == 0)
			throw runtime_error("No memory-saved simulation named " + mark);
		boost::iostreams::filtering_ostream out;
		if (boost::algorithm::ends_with(filename, ".bz2"))
			out.push(boost::iostreams::bzip2_compressor());
		out.push(boost::iostreams::file_sink(filename));
		if (!out.good())
			throw runtime_error("Error while opening file `" + filename + "' for writing.");
		LOG_INFO("Saving :memory:" << mark << " to " << filename);
		out << OMEGA.memSavedSimulations[":memory:" + mark];
	}
	string tmpToString(string mark)
	{
		if (OMEGA.memSavedSimulations.count(":memory:" + mark) == 0)
			throw runtime_error("No memory-saved simulation named " + mark);
		return OMEGA.memSavedSimulations[":memory:" + mark];
	}

	void reset()
	{
		OMEGA.stop();
		OMEGA.reset();
	}
	void resetThisScene()
	{
		Py_BEGIN_ALLOW_THREADS;
		OMEGA.stop();
		Py_END_ALLOW_THREADS;
		OMEGA.resetCurrentScene();
		OMEGA.createSimulationLoop();
	}
	void resetCurrentScene()
	{
		Py_BEGIN_ALLOW_THREADS;
		OMEGA.stop();
		Py_END_ALLOW_THREADS;
		OMEGA.resetCurrentScene();
		OMEGA.createSimulationLoop();
	}
	void resetTime()
	{
		OMEGA.getScene()->iter = 0;
		OMEGA.getScene()->time = 0;
		OMEGA.timeInit();
	}
	void switchScene() { std::swap(OMEGA.scenes[OMEGA.currentSceneNb], OMEGA.sceneAnother); }
	void resetAllScenes()
	{
		Py_BEGIN_ALLOW_THREADS;
		OMEGA.stop();
		Py_END_ALLOW_THREADS;
		OMEGA.resetAllScenes();
		OMEGA.createSimulationLoop();
	}
	void switchToScene(int i) { OMEGA.switchToScene(i); }
	int  addScene() { return OMEGA.addScene(); }

	// Scene manipulation in multithread situations as python object or as a string (e.g. FEMxDEM, MPI, ...)
	shared_ptr<Scene> scene_get() { return OMEGA.getScene(); }
	void              scene_set(const shared_ptr<Scene>& source)
	{
		Py_BEGIN_ALLOW_THREADS;
		reset();
		Py_END_ALLOW_THREADS;
		assertScene();
		OMEGA.setScene(source);
	}

#if PY_MAJOR_VERSION < 3
	string sceneToString() {
#else
	PyObject* sceneToString()
	{
#endif
		std::ostringstream oss;
	yade::ObjectIO::save<decltype(OMEGA.getScene()), boost::archive::binary_oarchive>(oss, "scene", OMEGA.getScene());
	oss.flush();
#if PY_MAJOR_VERSION < 3
	return oss.str();
#else
		const string s = oss.str();
		return PyBytes_FromStringAndSize(s.c_str(), s.length());
#endif
}

#ifdef YADE_MPI
//return vector<int> as contiguous bytes, skipping conversion to python list
PyObject*
intrsctToBytes(const shared_ptr<Subdomain>& subD, unsigned rank, bool mirror)
{
	if (subD->intersections.size() <= rank)
		LOG_ERROR("rank too large");
	if (not mirror)
		return PyBytes_FromStringAndSize((char*)&(subD->intersections[rank][0]), subD->intersections[rank].size() * sizeof(Body::id_t));
	else
		return PyBytes_FromStringAndSize((char*)&(subD->mirrorIntersections[rank][0]), subD->mirrorIntersections[rank].size() * sizeof(Body::id_t));
}
//return vector<int> as a writable contiguous array. Python syntax: bufferFromIntrsct(...)[:]=intrsctToBytes(...)
PyObject* bufferFromIntrsct(const shared_ptr<Subdomain>& subD, unsigned rank, unsigned size, bool mirror)
{
	//FIXME: if returning a memoryview, do we need to release it at some point? (https://docs.python.org/3/library/stdtypes.html#memoryview.release)
	if (subD->intersections.size() <= rank)
		LOG_ERROR("rank too large");
	if (not mirror) {
		if (subD->intersections[rank].size() != size)
			subD->intersections[rank].resize(size);
		return PyMemoryView_FromMemory((char*)(subD->intersections[rank].data()), subD->intersections[rank].size() * sizeof(Body::id_t), PyBUF_WRITE);
	} else {
		if (subD->mirrorIntersections[rank].size() != size)
			subD->mirrorIntersections[rank].resize(size);
		return PyMemoryView_FromMemory(
		        (char*)&(subD->mirrorIntersections[rank][0]), subD->mirrorIntersections[rank].size() * sizeof(Body::id_t), PyBUF_WRITE);
	}
}

#endif //YADE_MPI

void stringToScene(const string& sstring, string mark = "")
{
	Py_BEGIN_ALLOW_THREADS;
	OMEGA.stop();
	Py_END_ALLOW_THREADS;
	assertScene();
	OMEGA.memSavedSimulations[":memory:" + mark] = sstring;
	OMEGA.sceneFile                              = ":memory:" + mark;
	load(OMEGA.sceneFile, true);
}

int thisScene() { return OMEGA.currentSceneNb; }

void save(std::string fileName, bool quiet = false)
{
	assertScene();
	OMEGA.saveSimulation(fileName, quiet);
	// OMEGA.sceneFile=fileName; // done in Omega::saveSimulation;
}

py::list miscParams_get()
{
	py::list ret;
	FOREACH(shared_ptr<Serializable> & s, OMEGA.getScene()->miscParams) { ret.append(s); }
	return ret;
}

void miscParams_set(vector<shared_ptr<Serializable>> ss)
{
	vector<shared_ptr<Serializable>>& miscParams = OMEGA.getScene()->miscParams;
	miscParams.clear();
	FOREACH(shared_ptr<Serializable> s, ss) { miscParams.push_back(s); }
}


vector<shared_ptr<Engine>> engines_get(void)
{
	assertScene();
	Scene* scene = OMEGA.getScene().get();
	return scene->_nextEngines.empty() ? scene->engines : scene->_nextEngines;
}
void engines_set(const vector<shared_ptr<Engine>>& egs)
{
	assertScene();
	Scene* scene = OMEGA.getScene().get();
	if (scene->subStep < 0)
		scene->engines = egs; // not inside the engine loop right now, ok to update directly
	else
		scene->_nextEngines
		        = egs; // inside the engine loop, update _nextEngines; O.engines picks that up automatically, and Scene::moveToNextTimestep will put them in place of engines at the start of the next loop
	mapLabeledEntitiesToVariables();
}
// raw access to engines/_nextEngines, for debugging
vector<shared_ptr<Engine>> currEngines_get() { return OMEGA.getScene()->engines; }
vector<shared_ptr<Engine>> nextEngines_get() { return OMEGA.getScene()->_nextEngines; }

pyBodyContainer bodies_get(void)
{
	assertScene();
	return pyBodyContainer(OMEGA.getScene()->bodies);
}
pyInteractionContainer interactions_get(void)
{
	assertScene();
	return pyInteractionContainer(OMEGA.getScene()->interactions);
}

pyForceContainer    forces_get(void) { return pyForceContainer(OMEGA.getScene()); }
pyMaterialContainer materials_get(void) { return pyMaterialContainer(OMEGA.getScene()); }


py::list listChildClassesNonrecursive(const string& base)
{
	py::list ret;
	for (std::map<std::string, DynlibDescriptor>::const_iterator di = Omega::instance().getDynlibsDescriptor().begin();
	     di != Omega::instance().getDynlibsDescriptor().end();
	     ++di)
		if (Omega::instance().isInheritingFrom((*di).first, base))
			ret.append(di->first);
	return ret;
}

bool isChildClassOf(const string& child, const string& base) { return (Omega::instance().isInheritingFrom_recursive(child, base)); }

py::list plugins_get()
{
	const std::map<std::string, DynlibDescriptor>& plugins = Omega::instance().getDynlibsDescriptor();
	std::pair<string, DynlibDescriptor>            p;
	py::list                                       ret;
	for (const auto& p : plugins)
		ret.append(p.first);
	return ret;
}

pyTags tags_get(void)
{
	assertScene();
	return pyTags(OMEGA.getScene());
}

void interactionContainer_set(string clss)
{
	Scene* rb = OMEGA.getScene().get();
	if (rb->interactions->size() > 0)
		throw std::runtime_error("Interaction container not empty, will not change its class.");
	shared_ptr<InteractionContainer> ic = YADE_PTR_DYN_CAST<InteractionContainer>(ClassFactory::instance().createShared(clss));
	rb->interactions                    = ic;
}
string interactionContainer_get(string /*clss*/) { return OMEGA.getScene()->interactions->getClassName(); }

void bodyContainer_set(string clss)
{
	Scene* rb = OMEGA.getScene().get();
	if (rb->bodies->size() > 0)
		throw std::runtime_error("Body container not empty, will not change its class.");
	shared_ptr<BodyContainer> bc = YADE_PTR_DYN_CAST<BodyContainer>(ClassFactory::instance().createShared(clss));
	rb->bodies                   = bc;
}
string bodyContainer_get(string /*clss*/) { return OMEGA.getScene()->bodies->getClassName(); }
#ifdef YADE_OPENMP
int  numThreads_get() { return omp_get_max_threads(); }
void numThreads_set(int n)
{
	int bcn = OMEGA.getScene()->forces.getNumAllocatedThreads();
	if (bcn < n)
		LOG_WARN(
		        "ForceContainer has only " << bcn << " threads allocated. Changing thread number to on " << bcn << " instead of " << n
		                                   << " requested.");
	omp_set_num_threads(min(n, bcn));
	LOG_WARN("BUG: Omega().numThreads=n doesn't work as expected (number of threads is not changed globally). Set env var OMP_NUM_THREADS instead.");
}
#else
	int  numThreads_get() { return 1; }
	void numThreads_set(int) { LOG_WARN("Yade was compiled without openMP support, changing number of threads will have no effect."); }
#endif

shared_ptr<Cell> cell_get()
{
	if (OMEGA.getScene()->isPeriodic)
		return OMEGA.getScene()->cell;
	return shared_ptr<Cell>();
}
bool periodic_get(void) { return OMEGA.getScene()->isPeriodic; }
void periodic_set(bool v) { OMEGA.getScene()->isPeriodic = v; }

shared_ptr<EnergyTracker> energy_get() { return OMEGA.getScene()->energy; }
bool                      trackEnergy_get(void) { return OMEGA.getScene()->trackEnergy; }
void                      trackEnergy_set(bool e) { OMEGA.getScene()->trackEnergy = e; }

void disableGdb()
{
	signal(SIGSEGV, SIG_DFL);
	signal(SIGABRT, SIG_DFL);
}
void exitNoBacktrace(int status = 0)
{
	if (status == 0)
		signal(SIGSEGV, termHandlerNormal); /* unset the handler that runs gdb and prints backtrace */
	else
		signal(SIGSEGV, termHandlerError);
	// try to clean our mess
	Omega::instance().cleanupTemps();
	// flush all streams (so that in case we crash at exit, unflushed buffers are not lost)
	fflush(NULL);
	// attempt exit
	exit(status);
}
void runEngine(const shared_ptr<Engine>& e)
{
	LOG_WARN("Omega().runEngine(): deprecated, use __call__ method of the engine instance directly instead; will be removed in the future.");
	e->scene = OMEGA.getScene().get();
	e->action();
}
std::string tmpFilename() { return OMEGA.tmpFilename(); }
};

class pyGenericPotential : public GenericPotential, public py::wrapper<GenericPotential> {
public:
	Real potential(Real const& u, LubricationPhys const& phys) const
	{
		TRACE;
		return contactForce(u, phys.a) + potentialForce(u, phys.a);
	}

	void applyPotential(Real const& u, LubricationPhys& phys, Vector3r const& n)
	{
		TRACE;
		phys.normalContactForce   = contactForce(u, phys.a) * n;
		phys.normalPotentialForce = potentialForce(u, phys.a) * n;
		phys.contact              = hasContact(u, phys.a);
	}

	virtual Real contactForce(Real const& u, Real const& a) const
	{
		TRACE;
		gilLock lock;
#if PY_MAJOR_VERSION >= 3
		LOG_TRACE("GIL State: " << PyGILState_Check());
#endif
		return static_cast<Real>(get_override("contactForce")(u, a));
	}

	virtual Real potentialForce(Real const& u, Real const& a) const
	{
		TRACE;
		gilLock lock;
		return static_cast<Real>(get_override("potentialForce")(u, a));
	}

	virtual bool hasContact(Real const& u, Real const& a) const
	{
		TRACE;
		gilLock lock;
		return get_override("hasContact")(u, a);
	}
};

} // namespace yade


// BOOST_PYTHON_MODULE cannot be inside yade namespace, it has 'extern "C"' keyword, which strips it out of any namespaces.
BOOST_PYTHON_MODULE(wrapper)
try {
	using namespace yade; // 'using namespace' inside function keeps namespace pollution under control. Alernatively I could add y:: in front of function names below and put 'namespace y  = ::yade;' here.
	namespace py                = ::boost::python;
	py::scope().attr("__doc__") = "Wrapper for c++ internals of yade.";

	YADE_SET_DOCSTRING_OPTS;

	// keep in sync with lib/serialization/Serializable.hpp !
	py::enum_<yade::Attr::flags>("AttrFlags")
	        .value("noSave", yade::Attr::noSave)
	        .value("readonly", yade::Attr::readonly)
	        .value("triggerPostLoad", yade::Attr::triggerPostLoad)
	        .value("noResize", yade::Attr::noResize);

	py::class_<pyOmega>("Omega")
	        .add_property("iter", &pyOmega::iter, "Get current step number")
	        .add_property(
	                "subStep",
	                &pyOmega::subStep,
	                "Get the current subStep number (only meaningful if O.subStepping==True); -1 when outside the loop, otherwise either 0 "
	                "(O.subStepping==False) or number of engine to be run (O.subStepping==True)")
	        .add_property("subStepping", &pyOmega::subStepping_get, &pyOmega::subStepping_set, "Get/set whether subStepping is active.")
	        .add_property(
	                "stopAtIter", &pyOmega::stopAtIter_get, &pyOmega::stopAtIter_set, "Get/set number of iteration after which the simulation will stop.")
	        .add_property("stopAtTime", &pyOmega::stopAtTime_get, &pyOmega::stopAtTime_set, "Get/set time after which the simulation will stop.")
	        .add_property("time", &pyOmega::time, "Return virtual (model world) time of the simulation.")
	        .add_property("realtime", &pyOmega::realTime, "Return clock (human world) time the simulation has been running.")
	        .add_property("speed", &pyOmega::speed, "Return current calculation speed [iter/sec].")
	        .add_property("dt", &pyOmega::dt_get, &pyOmega::dt_set, "Current timestep (Δt) value.")
	        .add_property(
	                "dynDt",
	                &pyOmega::dynDt_get,
	                &pyOmega::dynDt_set,
	                "Whether a :yref:`TimeStepper` is used for dynamic Δt control. See :yref:`dt<Omega.dt>` on how to enable/disable :yref:`TimeStepper`.")
	        .add_property(
	                "dynDtAvailable",
	                &pyOmega::dynDtAvailable_get,
	                "Whether a :yref:`TimeStepper` is amongst :yref:`O.engines<Omega.engines>`, activated or not.")
	        .def("load",
	             &pyOmega::load,
	             (py::arg("file"), py::arg("quiet") = false),
	             "Load simulation from file. The file should be :yref:`saved<Omega.save>` in the same version of Yade, otherwise compatibility is not "
	             "guaranteed.")
	        .def("reload", &pyOmega::reload, (py::arg("quiet") = false), "Reload current simulation")
	        .def("save",
	             &pyOmega::save,
	             (py::arg("file"), py::arg("quiet") = false),
	             "Save current simulation to file (should be .xml or .xml.bz2 or .yade or .yade.gz). .xml files are bigger than .yade, but can be more or "
	             "less easily (due to their size) opened and edited, e.g. with text editors. .bz2 and .gz correspond both to compressed versions. All "
	             "saved files should be :yref:`loaded<Omega.load>` in the same version of Yade, otherwise compatibility is not guaranteed.")
	        .def("loadTmp",
	             &pyOmega::loadTmp,
	             (py::arg("mark") = "", py::arg("quiet") = false),
	             "Load simulation previously stored in memory by saveTmp. *mark* optionally distinguishes multiple saved simulations")
	        .def("saveTmp",
	             &pyOmega::saveTmp,
	             (py::arg("mark") = "", py::arg("quiet") = false),
	             "Save simulation to memory (disappears at shutdown), can be loaded later with loadTmp. *mark* optionally distinguishes different "
	             "memory-saved simulations.")
	        .def("lsTmp", &pyOmega::lsTmp, "Return list of all memory-saved simulations.")
	        .def("tmpToFile",
	             &pyOmega::tmpToFile,
	             (py::arg("fileName"), py::arg("mark") = ""),
	             "Save XML of :yref:`saveTmp<Omega.saveTmp>`'d simulation into *fileName*.")
	        .def("tmpToString", &pyOmega::tmpToString, (py::arg("mark") = ""), "Return XML of :yref:`saveTmp<Omega.saveTmp>`'d simulation as string.")
	        .def("run",
	             &pyOmega::run,
	             (py::arg("nSteps") = -1, py::arg("wait") = false),
	             "Run the simulation. *nSteps* how many steps to run, then stop (if positive); *wait* will cause not returning to python until simulation "
	             "will have stopped.")
	        .def("pause", &pyOmega::pause, "Stop simulation execution. (May be called from within the loop, and it will stop after the current step).")
	        .def("step", &pyOmega::step, "Advance the simulation by one step. Returns after the step will have finished.")
	        .def("wait", &pyOmega::wait, "Don't return until the simulation will have been paused. (Returns immediately if not running).")
	        .add_property("running", &pyOmega::isRunning, "Whether background thread is currently running a simulation.")
	        .add_property("filename", &pyOmega::get_filename, "Filename under which the current simulation was saved (None if never saved).")
	        .def("reset", &pyOmega::reset, "Reset simulations completely (including another scenes!).")
	        .def("resetThisScene", &pyOmega::resetThisScene, "Reset current scene.")
	        .def("resetCurrentScene", &pyOmega::resetCurrentScene, "Reset current scene.")
	        .def("resetAllScenes", &pyOmega::resetAllScenes, "Reset all scenes.")
	        .def("addScene", &pyOmega::addScene, "Add new scene to Omega, returns its number")
	        .def("switchToScene",
	             &pyOmega::switchToScene,
	             "Switch to defined scene. Default scene has number 0, other scenes have to be created by addScene method.")
	        .def("switchScene",
	             &pyOmega::switchScene,
	             "Switch to alternative simulation (while keeping the old one). Calling the function again switches back to the first one. Note that most "
	             "variables from the first simulation will still refer to the first simulation even after the switch\n(e.g. b=O.bodies[4]; "
	             "O.switchScene(); [b still refers to the body in the first simulation here])")
	        .add_property("thisScene", &pyOmega::thisScene, "Return current scene's id.")
	        .def("sceneToString",
	             &pyOmega::sceneToString,
	             "Return the entire scene as a string. Equivalent to using O.save(...) except that the scene goes to a string instead of a file. (see also "
	             "stringToScene())")
	        .def("stringToScene",
	             &pyOmega::stringToScene,
	             (py::arg("mark") = ""),
	             "Load simulation from a string passed as argument (see also sceneToString).")
	        .def("labeledEngine",
	             &pyOmega::labeled_engine_get,
	             // FIXME: use R"""(raw text)""" here, like exaplained in https://yade-dem.org/doc/prog.html#sphinx-documentation and used in py/_libVersions.cpp
	             "Return instance of engine/functor with the given label. This function shouldn't be called by the user directly; every ehange in "
	             "O.engines will assign respective global python variables according to labels.\n\nFor example:\n\n\t "
	             "*O.engines=[InsertionSortCollider(label='collider')]*\n\n\t *collider.nBins=5 # collider has become a variable after assignment to "
	             "O.engines automatically*")
	        .def("resetTime",
	             &pyOmega::resetTime,
	             "Reset simulation time: step number, virtual and real time. (Doesn't touch anything else, including timings).")
	        .def("plugins", &pyOmega::plugins_get, "Return list of all plugins registered in the class factory.")
	        // 		.def("_sceneObj",&pyOmega::scene_get,"Return the :yref:`scene <Scene>` object. Debugging only, all (or most) :yref:`Scene` functionality is proxies through :yref:`Omega`.")
	        .add_property(
	                "_sceneObj",
	                &pyOmega::scene_get,
	                &pyOmega::scene_set,
	                "Return the :yref:`scene <Scene>` object. Debugging only, all (or most) :yref:`Scene` functionality is proxies through :yref:`Omega`.")
	        .add_property(
	                "engines",
	                &pyOmega::engines_get,
	                &pyOmega::engines_set,
	                "List of engines in the simulation (corresponds to Scene::engines in C++ source code).")
	        .add_property("_currEngines", &pyOmega::currEngines_get, "Currently running engines; debugging only!")
	        .add_property(
	                "_nextEngines",
	                &pyOmega::nextEngines_get,
	                "Engines for the next step, if different from the current ones, otherwise empty; debugging only!")
	        .add_property(
	                "miscParams",
	                &pyOmega::miscParams_get,
	                &pyOmega::miscParams_set,
	                "MiscParams in the simulation (Scene::mistParams), usually used to save serializables that don't fit anywhere else, like GL functors")
	        .add_property("bodies", &pyOmega::bodies_get, "Bodies in the current simulation (container supporting index access by id and iteration)")
	        .add_property(
	                "interactions",
	                &pyOmega::interactions_get,
	                "Access to :yref:`interactions<Interaction>` of simulation, by using \n\n#. id's of both :yref:`Bodies<Body>` of the interactions, "
	                "e.g. ``O.interactions[23,65]``\n#. iteraction over the whole container::\n\n\tfor i in O.interactions: print i.id1,i.id2\n\n.. "
	                "note::\n\tIteration silently skips interactions that are not :yref:`real<Interaction.isReal>`.")
	        .add_property("materials", &pyOmega::materials_get, "Shared materials; they can be accessed by id or by label")
	        .add_property("forces", &pyOmega::forces_get, ":yref:`ForceContainer` (forces, torques, displacements) in the current simulation.")
	        .add_property(
	                "energy",
	                &pyOmega::energy_get,
	                ":yref:`EnergyTracker` of the current simulation. (meaningful only with :yref:`O.trackEnergy<Omega.trackEnergy>`)")
	        .add_property(
	                "trackEnergy", &pyOmega::trackEnergy_get, &pyOmega::trackEnergy_set, "When energy tracking is enabled or disabled in this simulation.")
	        .add_property(
	                "tags",
	                &pyOmega::tags_get,
	                "Tags (string=string dictionary) of the current simulation (container supporting string-index access/assignment)")
	        .def("childClassesNonrecursive",
	             &pyOmega::listChildClassesNonrecursive,
	             "Return list of all classes deriving from given class, as registered in the class factory")
	        .def("isChildClassOf", &pyOmega::isChildClassOf, "Tells whether the first class derives from the second one (both given as strings).")
	        .add_property(
	                "timingEnabled",
	                &pyOmega::timingEnabled_get,
	                &pyOmega::timingEnabled_set,
	                "Globally enable/disable timing services (see documentation of the :yref:`timing module<yade.timing>`).")
	        .add_property(
	                "forceSyncCount",
	                &pyOmega::forceSyncCount_get,
	                &pyOmega::forceSyncCount_set,
	                "Counter for number of syncs in ForceContainer, for profiling purposes.")
	        .add_property("numThreads", &pyOmega::numThreads_get /* ,&pyOmega::numThreads_set*/, "Get maximum number of threads openMP can use.")
	        .add_property("cell", &pyOmega::cell_get, "Periodic cell of the current scene (None if the scene is aperiodic).")
	        .add_property("periodic", &pyOmega::periodic_get, &pyOmega::periodic_set, "Get/set whether the scene is periodic or not (True/False).")
	        .def("exitNoBacktrace",
	             &pyOmega::exitNoBacktrace,
	             (py::arg("status") = 0),
	             "Disable SEGV handler and exit, optionally with given status number.")
	        .def("disableGdb", &pyOmega::disableGdb, "Revert SEGV and ABRT handlers to system defaults.")
	        .def("runEngine",
	             &pyOmega::runEngine,
	             "Run given engine exactly once; simulation time, step number etc. will not be incremented (use only if you know what you do).")
	        .def("tmpFilename", &pyOmega::tmpFilename, "Return unique name of file in temporary directory which will be deleted when yade exits.")
#ifdef YADE_MPI
	        .def("intrsctToBytes",
	             &pyOmega::intrsctToBytes,
	             (py::arg("subdomain"), py::arg("rank"), py::arg("mirror")),
	             "returns a copy of intersections[rank] (a vector<int>) from a subdomain in the form of bytes. Returns a copy mirrorIntersections[rank] if "
	             "mirror=True.")
	        .def("bufferFromIntrsct",
	             &pyOmega::bufferFromIntrsct,
	             (py::arg("subdomain"), py::arg("rank"), py::arg("size"), py::arg("mirror")),
	             "returns a (char*) pointer to the underying buffer of intersections[rank], so that it can be overwritten. Size must be passed in advance. "
	             "Pointer to mirrorIntersections[rank] is returned if mirror=True. Python syntax: bufferFromIntrsct(...)[:]=bytes(something)")
#endif
	        ;
	py::class_<pyTags>(
	        "TagsWrapper",
	        "Container emulating dictionary semantics for accessing tags associated with simulation. Tags are accesed by strings.",
	        py::init<pyTags&>())
	        .def("__getitem__", &pyTags::getItem)
	        .def("__setitem__", &pyTags::setItem)
	        .def("keys", &pyTags::keys)
	        .def("has_key", &pyTags::hasKey);
	py::class_<pyBodyContainer>("BodyContainer", py::init<pyBodyContainer&>())
	        .def("__getitem__", &pyBodyContainer::pyGetitem)
	        .def("__len__", &pyBodyContainer::length)
	        .def("__iter__", &pyBodyContainer::pyIter)
	        .def("append", &pyBodyContainer::append, "Append one Body instance, return its id.")
	        .def("append", &pyBodyContainer::appendList, "Append list of Body instance, return list of ids")
	        .def("appendClumped",
	             &pyBodyContainer::appendClump,
	             (py::arg("discretization") = 0),
	             "Append given list of bodies as a clump (rigid aggregate); returns a tuple of ``(clumpId,[memberId1,memberId2,...])``. Clump masses and "
	             "inertia are adapted automatically (for details see :yref:`clump()<BodyContainer.clump>`).")
	        .def("clump",
	             &pyBodyContainer::clump,
	             (py::arg("discretization") = 0),
	             "Clump given bodies together (creating a rigid aggregate); returns ``clumpId``. Clump masses and inertia are adapted automatically when "
	             "discretization>0. If clump members are overlapping this is done by integration/summation over mass points using a regular grid of cells "
	             "(grid cells length is defined as $R_{min}/discretization$, where $R_{min}$ is minimum clump member radius). For non-overlapping members "
	             "inertia of the clump is the sum of inertias from members. If discretization<=0 sum of inertias from members is used (faster, but "
	             "inaccurate).")
	        .def("updateClumpProperties",
	             &pyBodyContainer::updateClumpProperties,
	             (py::arg("excludeList") = py::list(), py::arg("discretization") = 5),
	             "Manually force Yade to update clump properties mass, volume and inertia (for details of 'discretization' value see "
	             ":yref:`clump()<BodyContainer.clump>`). Can be used, when clumps are modified or erased during a simulation. Clumps can be excluded from "
	             "the calculation by giving a list of ids: *O.bodies.updateProperties([ids])*.")
	        .def("deleteClumpMember", &pyBodyContainer::deleteClumpMember, "Erase clump member.") //FIXME
	        .def("deleteClumpBody", &pyBodyContainer::deleteClumpBody, "Erase clump member.")     //FIXME
	        .def("addToClump",
	             &pyBodyContainer::addToClump,
	             (py::arg("discretization") = 0),
	             "Add body b (or a list of bodies) to an existing clump c. c must be clump and b may not be a clump member of c. Clump masses and inertia "
	             "are adapted automatically (for details see :yref:`clump()<BodyContainer.clump>`).\n\nSee :ysrc:`examples/clumps/addToClump-example.py` "
	             "for an example script.\n\n.. note:: If b is a clump itself, then all members will be added to c and b will be deleted. If b is a clump "
	             "member of clump d, then all members from d will be added to c and d will be deleted. If you need to add just clump member b, "
	             ":yref:`release<BodyContainer.releaseFromClump>` this member from d first.")
	        .def("releaseFromClump",
	             &pyBodyContainer::releaseFromClump,
	             (py::arg("discretization") = 0),
	             "Release body b from clump c. b must be a clump member of c. Clump masses and inertia are adapted automatically (for details see "
	             ":yref:`clump()<BodyContainer.clump>`).\n\nSee :ysrc:`examples/clumps/releaseFromClump-example.py` for an example script.\n\n.. note:: If "
	             "c contains only 2 members b will not be released and a warning will appear. In this case clump c should be "
	             ":yref:`erased<BodyContainer.erase>`.")
	        .def("replaceByClumps",
	             &pyBodyContainer::replaceByClumps,
	             (py::arg("discretization") = 0),
	             "Replace spheres by clumps using a list of clump templates and a list of amounts; returns a list of tuples: "
	             "``[(clumpId1,[memberId1,memberId2,...]),(clumpId2,[memberId1,memberId2,...]),...]``. A new clump will have the same volume as the "
	             "sphere, that was replaced. Clump masses and inertia are adapted automatically (for details see :yref:`clump()<BodyContainer.clump>`). "
	             "\n\n\t *O.bodies.replaceByClumps( [utils.clumpTemplate([1,1],[.5,.5])] , [.9] ) #will replace 90 % of all standalone spheres by "
	             "'dyads'*\n\nSee :ysrc:`examples/clumps/replaceByClumps-example.py` for an example script.")
	        .def("getRoundness",
	             &pyBodyContainer::getRoundness,
	             (py::arg("excludeList") = py::list()),
	             "Returns roundness coefficient RC = R2/R1. R1 is the equivalent sphere radius of a clump. R2 is the minimum radius of a sphere, that "
	             "imbeds the clump. If just spheres are present RC = 1. If clumps are present 0 < RC < 1. Bodies can be excluded from the calculation by "
	             "giving a list of ids: *O.bodies.getRoundness([ids])*.\n\nSee :ysrc:`examples/clumps/replaceByClumps-example.py` for an example script.")
	        .def("clear", &pyBodyContainer::clear, "Remove all bodies (interactions not checked)")
	        .def("erase",
	             &pyBodyContainer::erase,
	             (py::arg("eraseClumpMembers") = 0),
	             "Erase body with the given id; all interaction will be deleted by InteractionLoop in the next step. If a clump is erased use "
	             "*O.bodies.erase(clumpId,True)* to erase the clump AND its members.")
	        .def("replace", &pyBodyContainer::replace)
	        .def("insertAtId", &pyBodyContainer::insertAtId, (py::arg("insertatid")), "Insert a body at theid, (no body should exist in this id)")
#ifdef YADE_MPI
	        .def("subdomainBodies", &pyBodyContainer::subdomainBodies, "id's of bodies with bounds in MPI subdomain")
#endif //YADE_MPI
	        .add_property(
	                "useRedirection",
	                &pyBodyContainer::getUseRedirection,
	                &pyBodyContainer::setUseRedirection,
	                "true if the scene uses up-to-date lists for boundedBodies and realBodies; turned true automatically 1/ after removal of bodies if "
	                ":yref:`enableRedirection`=True, and 2/ in MPI execution.")
	        .add_property(
	                "enableRedirection",
	                &pyBodyContainer::getEnableRedirection,
	                &pyBodyContainer::setEnableRedirection,
	                "let collider switch to optimized algorithm with body redirection when bodies are erased - true by default");
	py::class_<pyBodyIterator>("BodyIterator", py::init<pyBodyIterator&>())
	        .def("__iter__", &pyBodyIterator::pyIter)
	        .def("__next__", &pyBodyIterator::pyNext);
	py::class_<pyInteractionContainer>(
	        "InteractionContainer",
	        "Access to :yref:`interactions<Interaction>` of simulation, by using \n\n#. id's of both :yref:`Bodies<Body>` of the interactions, e.g. "
	        "``O.interactions[23,65]``\n#. iteraction over the whole container::\n\n\tfor i in O.interactions: print i.id1,i.id2\n\n.. note::\n\tIteration "
	        "silently skips interactions that are not :yref:`real<Interaction.isReal>`.",
	        py::init<pyInteractionContainer&>())
	        .def("__iter__", &pyInteractionContainer::pyIter)
	        .def("__getitem__", &pyInteractionContainer::pyGetitem)
	        .def("__len__", &pyInteractionContainer::len)
	        .def("has", &pyInteractionContainer::has, "Tell if a pair of ids corresponds to an existing interaction (real or not)")
	        .def("countReal", &pyInteractionContainer::countReal, "Return number of interactions that are \"real\", i.e. they have phys and geom.")
	        .def("nth",
	             &pyInteractionContainer::pyNth,
	             "Return n-th interaction from the container (usable for picking random interaction). The virtual interactions are not reached.")
	        .def("withBody", &pyInteractionContainer::withBody, "Return list of real interactions of given body.")
	        .def("withBodyAll", &pyInteractionContainer::withBodyAll, "Return list of all (real as well as non-real) interactions of given body.")
	        .def("all",
	             &pyInteractionContainer::getAll,
	             (py::arg("onlyReal") = false),
	             "Return list of all interactions. Virtual interaction are filtered out if onlyReal=True, else (default) it dumps the full content.")
	        .def("eraseNonReal", &pyInteractionContainer::eraseNonReal, "Erase all interactions that are not :yref:`real <Interaction.isReal>`.")
	        .def("erase",
	             &pyInteractionContainer::erase,
	             "Erase one interaction, given by id1, id2 (internally, ``requestErase`` is called -- the interaction might still exist as potential, if "
	             "the :yref:`Collider` decides so).")
	        .add_property("serializeSorted", &pyInteractionContainer::serializeSorted_get, &pyInteractionContainer::serializeSorted_set)
	        .def("clear",
	             &pyInteractionContainer::clear,
	             "Remove all interactions, and invalidate persistent collider data (if the collider supports it).");
	py::class_<pyInteractionIterator>("InteractionIterator", py::init<pyInteractionIterator&>())
	        .def("__iter__", &pyInteractionIterator::pyIter)
	        .def("__next__", &pyInteractionIterator::pyNext);

	py::class_<pyForceContainer>("ForceContainer", py::init<pyForceContainer&>())
	        .def("f",
	             &pyForceContainer::force_get,
	             (py::arg("id"), py::arg("sync") = false),
	             "Resultant force on body, excluding :yref:`gravity<NewtonIntegrator.gravity>`. For clumps in openMP, synchronize the force container with "
	             "sync=True, else the value will be wrong.")
	        .def("t",
	             &pyForceContainer::torque_get,
	             (py::arg("id"), py::arg("sync") = false),
	             "Torque applied on body. For clumps in openMP, synchronize the force container with sync=True, else the value will be wrong.")
	        .def("m", &pyForceContainer::torque_get, (py::arg("id"), py::arg("sync") = false), "Deprecated alias for t (torque).")
	        .def("addF",
	             &pyForceContainer::force_add,
	             (py::arg("id"), py::arg("f"), py::arg("permanent") = false),
	             "Apply force on body (accumulates). The force applies for one iteration, then it is reset by ForceResetter. \n # permanent parameter is "
	             "deprecated, instead of addF(...,permanent=True) use setPermF(...).")
	        .def("addT",
	             &pyForceContainer::torque_add,
	             (py::arg("id"), py::arg("t"), py::arg("permanent") = false),
	             "Apply torque on body (accumulates). The torque applies for one iteration, then it is reset by ForceResetter. \n # permanent parameter is "
	             "deprecated, instead of addT(...,permanent=True) use setPermT(...).")
	        .def("permF", &pyForceContainer::permForce_get, (py::arg("id")), "read the value of permanent force on body (set with setPermF()).")
	        .def("permT", &pyForceContainer::permTorque_get, (py::arg("id")), "read the value of permanent torque on body (set with setPermT()).")
	        .def("setPermF", &pyForceContainer::permForce_set, "set the value of permanent force on body.")
	        .def("setPermT", &pyForceContainer::permTorque_set, "set the value of permanent torque on body.")
	        .def("reset",
	             &pyForceContainer::reset,
	             (py::arg("resetAll") = true),
	             "Reset the force container, including user defined permanent forces/torques. resetAll=False will keep permanent forces/torques unchanged.")
	        .def("getPermForceUsed", &pyForceContainer::getPermForceUsed, "Check wether permanent forces are present.")
	        .add_property(
	                "syncCount",
	                &pyForceContainer::syncCount_get,
	                &pyForceContainer::syncCount_set,
	                "Number of synchronizations  of ForceContainer (cummulative); if significantly higher than number of steps, there might be unnecessary "
	                "syncs hurting performance.");

	py::class_<pyMaterialContainer>(
	        "MaterialContainer",
	        "Container for :yref:`Materials<Material>`. A material can be accessed using \n\n #. numerical index in range(0,len(cont)), like cont[2]; \n "
	        "#. textual label that was given to the material, like cont['steel']. This entails traversing all materials and should not be used frequently.",
	        py::init<pyMaterialContainer&>())
	        .def("append", &pyMaterialContainer::append, "Add new shared :yref:`Material`; changes its id and return it.")
	        .def("append", &pyMaterialContainer::appendList, "Append list of :yref:`Material` instances, return list of ids.")
	        .def("index", &pyMaterialContainer::index, "Return id of material, given its label.")
	        .def("__getitem__", &pyMaterialContainer::getitem_id)
	        .def("__getitem__", &pyMaterialContainer::getitem_label)
	        .def("__len__", &pyMaterialContainer::len);

	py::class_<STLImporter>("STLImporter").def("ymport", &STLImporter::import);

	py::class_<pyGenericPotential, boost::noncopyable>("GenericPotential")
	        .def("contactForce", py::pure_virtual(&pyGenericPotential::contactForce), "This function should return contact force norm.")
	        .def("potentialForce", py::pure_virtual(&pyGenericPotential::potentialForce), "This function should return potential force norm.")
	        .def("hasContact", py::pure_virtual(&pyGenericPotential::hasContact), "This function should return true if there are contact.");

	//////////////////////////////////////////////////////////////
	///////////// proxyless wrappers
	Serializable().pyRegisterClass(py::scope());

	py::class_<TimingDeltas, shared_ptr<TimingDeltas>, boost::noncopyable>("TimingDeltas")
	        .add_property("data", &TimingDeltas::pyData, "Get timing data as list of tuples (label, execTime[nsec], execCount) (one tuple per checkpoint)")
	        .def("reset", &TimingDeltas::reset, "Reset timing information");

	py::scope().attr("O") = pyOmega();

} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}
