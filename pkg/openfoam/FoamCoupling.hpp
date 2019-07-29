// YADE-OpenFOAM coupling module, Deepak kn  deepak.kunhappan@3sr-grenoble.fr/deepak.kn1990@gmail.com
#ifdef YADE_MPI

#pragma once 

#include <core/Scene.hpp> 
#include <core/GlobalEngine.hpp> 
#include <core/Body.hpp>
#include <core/Omega.hpp>
#include <pkg/common/Sphere.hpp> 
#include <vector> 
#include <core/InteractionContainer.hpp> 
#include <lib/serialization/Serializable.hpp>



namespace yade { // Cannot have #include directive inside.

#include <core/Subdomain.hpp>
#include <pkg/common/Aabb.hpp>
#include <pkg/common/Dispatching.hpp>
class Scene; 

class FoamCoupling : public GlobalEngine {
	private:

	// some variables for MPI_Send/Recv 
		const int sendTag=500;  
		MPI_Status status; 
		int rank, commSize; 
		int szdff, localCommSize, worldCommSize; 

  private:
    // some variables for MPI_Send/Recv 
    const int sendTag=500;  
    MPI_Status status; 
    int rank, commSize; 


  public: 
    
    void getRank(); 
    void setNumParticles(int); 
    void setIdList(const std::vector<int>& );  
    void killMPI(); 
    void updateProcList();
    void castParticle();
    void castNumParticle(int); 
    void castTerminate();  
    void resetProcList(); 
    void recvHydroForce(); 
    void setHydroForce();
    void sumHydroForce(); 
    void exchangeDeltaT();  
    void runCoupling(); 
    bool exchangeData();
    Real getViscousTimeScale();  // not fully implemented, piece of code left in foam.
    virtual void action(); 
    virtual ~FoamCoupling(){}; 
    std::vector<int> bodyList; 
    std::vector<double> hydroForce; 
    std::vector<double> particleData;
    std::vector<int>  procList; 
    Real foamDeltaT; 
    long int  dataExchangeInterval=1; 
    bool recvdFoamDeltaT; 
    bool isGaussianInterp;
    bool initDone = false; 
    void insertBodyId(int); 
    bool eraseId(int);
    int getNumBodies(); 
    std::vector<int> getIdList(); 
    
		void getRank(); 
		void setNumParticles(int); 
		void setIdList(const boost::python::list& );  
		void killMPI(); 
		void updateProcList();
		void castParticle();
		void castNumParticle(int); 
		void resetProcList(); 
		void recvHydroForce(); 
		void setHydroForce();
		void sumHydroForce(); 
		void exchangeDeltaT();  
		void runCoupling(); 
		bool exchangeData();
                void castTerminate(); 
		Real getViscousTimeScale();  // not fully implemented, piece of code left in foam.
		virtual void action(); 
		virtual ~FoamCoupling(){}; 
		std::vector<int> bodyList; 
		std::vector<double> hydroForce; 
		std::vector<double> particleData;
		std::vector<int>  procList; 
		std::vector<Body::id_t> fluidDomains; 
		std::vector<std::vector<Body::id_t> > inBoxIds; 
		Real foamDeltaT; 
		long int  dataExchangeInterval=1; 
		bool recvdFoamDeltaT; 
		bool isGaussianInterp;
		void getFluidDomainBbox(); 
		
      
    YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(FoamCoupling,GlobalEngine, "An engine for coupling Yade with the finite volume fluid solver OpenFOAM in parallel." " \n Requirements : Yade compiled with MPI libs, OpenFOAM-6 (openfoam is not required for compilation)." "Yade is executed under MPI environment with OpenFOAM simultaneously, and using MPI communication  routines data is exchanged between the solvers."
   " \n \n 1. Yade broadcasts the particle data -> position, velocity, ang-velocity, radius to all the foam processes as in :yref:`castParticle <FoamCoupling::castParticle>` \n"
  "2. In each foam process, particle is searched.Yade keeps a vector(cpp) of the rank of the openfoam process containing that particular particle (FoamCoupling::procList), using :yref:`updateProcList <FoamCoupling::updateProcList>`\n"
  "3. In simple lagrangian point force coupling Yade recieves the particle hydrodynamic force and torque from the openfoam process, the sender is identified from  the vector :yref:`FoamCoupling::procList`.\n"
  "In the case of Gaussian interpolation, contribution from every  process is summed using function :yref:`sumHydroForce <FoamCoupling::sumHydroForce>`. \n"
  "4. The interval (substepping) is set automatically (:yref:`FoamCoupling::dataExchangeInterval`) based on dtfoam/dtYade, calculated in function :yref:`exchangeDeltaT<FoamCoupling::exchangeDeltaT>`  ", 
    ((int,numParticles,1, , "number of particles in coupling."))
    ((double,particleDensity,1, , "particle Density")) //not needed  as this is set in foam  
    ((double,fluidDensity,1, ,"fluidDensity")) //not needed  as this is set in foam  
    ,
    ,
    ,
    .def("setIdList", &FoamCoupling::setIdList,boost::python::arg("bodyIdlist"), "list of body ids in hydroForce coupling. (links to :yref: `FoamCoupling::bodyList` vector, used to build particle data :yref:`FoamCoupling::particleData`. :yref:`FoamCoupling::particleData` contains the particle pos, vel, angvel, radius and this is sent to foam. )")
    .def("getRank", &FoamCoupling::getRank, "Initiallize MPI communicator for coupling. Should be called at the beginning of the script. :yref: `initMPI <FoamCoupling::initMPI>` Initializes  the MPI environment. " )
    .def("killMPI", &FoamCoupling::killMPI, "Destroy MPI, to be called at the end of the simulation, from :yref:`killMPI<FoamCoupling::killMPI>`") 
    .def("setNumParticles",&FoamCoupling::setNumParticles,boost::python::arg("numparticles"),"number of particles in coupling")
    .def("insertBodyId", &FoamCoupling::insertBodyId,boost::python::arg("newId"), "insert a new body id for hydrodynamic force coupling") 
    .def("eraseId", &FoamCoupling::eraseId,boost::python::arg("idToErase"), "remove a body from hydrodynamic force coupling")
    .def("getNumBodies", &FoamCoupling::getNumBodies, "get the number of bodies in the coupling") 
    .def("getIdList", &FoamCoupling::getIdList, "get the ids of bodies in coupling")
    .def_readonly("foamDeltaT", &FoamCoupling::foamDeltaT, "timestep in openfoam solver from  :yref:`exchangeDeltaT <FoamCoupling::exchangeDeltaT>` ") 
    .def_readonly("dataExchangeInterval", &FoamCoupling::dataExchangeInterval, "Number of iterations/substepping : for stability and to be in sync with fluid solver calculated in :yref:`exchangeDeltaT <FoamCoupling::exchangeDeltaT>`")
    .def_readwrite("isGaussianInterp", &FoamCoupling::isGaussianInterp, "switch for Gaussian interpolation of field varibles in openfoam. Uses  :yref:`sumHydroForce<FoamCoupling::sumHydroForce>` to obtain hydrodynamic force ") 
    )
    DECLARE_LOGGER; 
}; 
REGISTER_SERIALIZABLE(FoamCoupling); 



/* a class for holding info on the min and max bounds of the fluid mesh. Each fluid proc has a domain minmax and */
class FluidDomainBbox : public Shape{
	public:
		std::vector<double> minMaxBuff; // a buffer to receive the min max during MPI_Recv. 
		void setMinMax(std::vector<double>& minmaxbuff){
			if (minMaxBuff.size() != 6){LOG_ERROR("incorrect minmaxbuff size. FAIL"); return; }
			minBound[0] = minmaxbuff[0]; 
			minBound[1] = minmaxbuff[1];
			minBound[2] = minmaxbuff[2]; 
			maxBound[0] = minmaxbuff[3]; 
			maxBound[1] = minmaxbuff[4]; 
			maxBound[2] = minmaxbuff[5]; 
			minMaxisSet = true; 
			
		}
		virtual ~FluidDomainBbox() {}; 
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(FluidDomainBbox,Shape,"The bounding box of a fluid grid from one OpenFOAM/YALES2 proc",
		((int,domainRank,-1,,"rank of the OpenFOAM/YALES2 proc"))
		((bool,minMaxisSet,false,,"flag to check if the min max bounds of this body are set."))
		((std::vector<Body::id_t>,bIds,std::vector<Body::id_t>(), , "ids of bodies intersecting with this subdomain, "))
		((Vector3r,minBound,Vector3r(NaN, NaN, NaN),,"min bounds of the fluid grid ")) 
		((Vector3r,maxBound,Vector3r(NaN, NaN, NaN),,"max bounds of the fluid grid"))	
		((bool,hasIntersection,false,,"if this Yade subdomain has intersection with this OpenFOAM/YALES2 subdomain"))
		,
		createIndex(); 
		,
		); 
		DECLARE_LOGGER;
		REGISTER_CLASS_INDEX(FluidDomainBbox, Shape);   
}; 
REGISTER_SERIALIZABLE(FluidDomainBbox); 

class Bo1_FluidDomainBbox_Aabb : public BoundFunctor{
	public:
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*); 
		FUNCTOR1D(FluidDomainBbox);  
		YADE_CLASS_BASE_DOC(Bo1_FluidDomainBbox_Aabb,BoundFunctor, "creates/updates an :yref:`Aabb` of a :yref:`FluidDomainBbox`."); 
}; 
REGISTER_SERIALIZABLE(Bo1_FluidDomainBbox_Aabb); 



#endif  
}
