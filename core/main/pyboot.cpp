#include <lib/base/Logging.hpp>
#include <core/Omega.hpp>

#include <signal.h>

#ifdef YADE_DEBUG
void crashHandler(int sig)
{
	switch (sig) {
		case SIGABRT:
		case SIGSEGV:
			signal(SIGSEGV, SIG_DFL);
			signal(SIGABRT, SIG_DFL); // prevent loops - default handlers
			std::cerr << "SIGSEGV/SIGABRT handler called; gdb batch file is `" << yade::Omega::instance().gdbCrashBatch << "'" << std::endl;
			int ret = std::system((std::string("gdb -x ") + yade::Omega::instance().gdbCrashBatch).c_str());
			if (ret != 0)
				std::cerr << "unable to execute gdb" << std::endl;
			raise(sig); // reemit signal after exiting gdb
			break;
	}
}
#endif

namespace forCtags {
struct pyboot {
}; // for ctags
}

/* Initialize yade, loading given plugins */
void yadeInitialize(boost::python::list& pp, const std::string& confDir)
{
	PyEval_InitThreads();

	yade::Omega& O(yade::Omega::instance());
	O.init();
	O.origArgv = NULL;
	O.origArgc = 0; // not needed, anyway
	O.confDir  = confDir;
	O.initTemps();
#ifdef YADE_DEBUG
	std::ofstream gdbBatch;
	O.gdbCrashBatch = O.tmpFilename();
	gdbBatch.open(O.gdbCrashBatch.c_str());
	gdbBatch << "attach " << boost::lexical_cast<std::string>(getpid()) << "\nset pagination off\nthread apply all backtrace\ndetach\nquit\n";
	gdbBatch.close();
	signal(SIGABRT, crashHandler);
	signal(SIGSEGV, crashHandler);
#endif

	std::vector<std::string> ppp;
	for (int i = 0; i < boost::python::len(pp); i++) {
		ppp.push_back(boost::python::extract<std::string>(pp[i]));
	}
	yade::Omega::instance().loadPlugins(ppp);
}
void yadeFinalize() { yade::Omega::instance().cleanupTemps(); }

BOOST_PYTHON_MODULE(boot)
{
	boost::python::scope().attr("initialize") = &yadeInitialize;
	boost::python::scope().attr("finalize")   = &yadeFinalize; //,"Finalize yade (only to be used internally).")
}
