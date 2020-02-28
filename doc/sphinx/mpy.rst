.. _mpy:

MPI parallelization
===================

The module :yref:`yade.mpy` implements parallelization by domain decomposition (distributed memory) using the Message Passing Interface (MPI) implemented by OpenMPI. It aims at exploiting large numbers of cores, where shared memory techniques (OpenMP) are helpess.
The shared memory and the distributed memory approaches are compatible, i.e. it is possible to run hybrid jobs using both, as shown below (and it may well be the optimal solution in many cases).

Most calls to OpenMPI library are done in Python using `mpi4py <https://mpi4py.readthedocs.io>`_. For the sake of efficiency some critical communications are triggered via python wrappers of C++ functions, wherein messages are produced, sent/received, and processed.

.. note:: Disclaimer: even though the `yade.mpy` module provides the function :yref:`mpirun<yade.mpy.mpirun>`, which may seem as a simple replacement for `O.run()`, setting up a simulation with mpy might be deceptively triavial.
    As of now, it is anticipated that, in general, a simple replacement of "run" by "mpirun" in an arbitrary script will not speedup anything and may even fail miserably (it could be improved in the future). To understand why, and to tackle the problems, basic knowledge of how MPI works will certainly help (specifically `mpi4py <https://mpi4py.readthedocs.io>`_), and carefull reading of this page is mandatory.
    Suggestions on how to improve this documentation or the implementation are welcome. Uninformed questions might be simply ignored and be refered to this page for the time being.


Concepts
________

**subdomain**: a (sub)set of bodies attached to one MPI process after domain decomposition. The corresponding class in Yade is :yref:`Subdomain`, a `Shape` instance with helper functions for MPI communications. In some sense `Subdomain` is to subscribed bodies what :yref:`Clump` (another `Shape`) is to clump members.

**rank**: subdomain index from 0 to *N*-1  (with *N* the number of mpi processes) to identify subdomains. The rank of the subdomain a body belongs to can be retrieved as :yref:`Body.subdomain`. Each subdomain corresponds to an instance of yade and a specific scene during parallel execution. The rank of the scene is given by :yref:`Scene.subdomain`.

**master**: refers to subdomain with *rank* =0. This subdomain does not behave like others. In general master will handle boundary conditions and it will control transitions and termination of the whole simulation. Unlike standard subdomains it may not contain a large number of raw bodies (i.e. not beyond objects bounding the scene such as walls or boxes). In interactive execution master is the process responding to the python prompt.

**intersections**: subsets of bodies in a subdomain intersected by the bounding box of other subdomains (see `fig-subdomains`_). *intersection(i,j)* refers to the bodies owned by current (*i*) subdomain and intersecting subdomain *j* (retrieved as :yref:`O._sceneObj.subD.intersections[j]<Subdomain.intersections>`); *mirrorIntersection(i,j)* refers to bodies owned by *j* and intersecting current domain (retrieved as :yref:`O._sceneObj.subD.mirrorIntersections[j]<Subdomain.mirrorIntersections>`). The bodies are listed by :yref:`Body.id`. By definition *intersection(i,j)=mirrorIntersection(j,i)*.

The intersections and mirror intersections are updated automatically as part of parallel collision detection. They define which body states need to be communicated. The bodies in intersections need to be *sent* to other subdomains (in pratice only updated position and velocity are sent at every iteration), the bodies in mirrorIntersections need to be received from other subdomains.


.. _fig-subdomains:
.. figure:: fig/subdomains.*
	:width: 12cm

Two overlapping subdomains and their intersections. In this situation we have *SubD1.intersections[SubD2.subdomain]=[id4,id5]* and *SubD1.mirrorIntersections[SubD2.subdomain]=[id1]*, with *SubD1* and *SubD2* instances of :yref:`Subdomain`.


.. _sect_mpi_implementation:

Implementation
---------------

For demonstrating the main internal steps in the implemented parallel algorithm let us conider the example script :ysrc:`examples/mpi/testMPI_2D.py`. Executing this script (interactive or passive mode) with three MPI processes generates the scene as shown in `fig-scene`_. It then executes :yref:`mpirun<yade.mpy.mpirun>`, which triggers the steps described hereafter.

.. _fig-scene:
.. figure:: fig/mpy_schema0.*
	:width: 25%
	:align: center


In this scene, we have three MPI processes (three subdomains) and the raw bodies are partitioned among the subdomains/ranks 1 and 2. The master process with subdomain=0 holds the boundary/wall type body. Bodies can be manually assigned or automatically assigned via a domain decomposition algorithm. Details 
on the dommain decomposition algorithm is presented in the later section of this document. 

**Scene splitting** :

In the function :yref:`yade.mpy.splitScene`, called at the beginning of mpi execution, specific engines are added silently to the scene in order to handle what will happen next. That very intrusive operation can even change settings of some pre-existing engines, in particular :yref:`InsertionSortCollider`, to make them behave with MPI-friendlyness. :yref:`InsertionSortCollider.verletDist` is an important factor controlling the efficiency of the simulations. The reason for this will become evident in the later steps. 

**Bounds dispatching** : In the next step, the :yref:`Body.bound` is dispatched with the :yref:`Aabb` extended as shown in figure `fig-regularbounds`_ (in dotted lines). Note that the :yref:`Subdomain` :yref:`Aabb` is obtained from taking the min and max of the owned bodies, see figure `fig-subDBounds`_  
with solid coloured lines for the subdomain :yref:`Aabb`. At this time, the min and max of other subdomains are unknown. 

.. _fig-regularbounds:
.. figure:: fig/mpy_schema1a.*
	:width: 25%
	:align: center


.. _fig-subDBounds:
.. figure:: fig/mpy_schema1b.*
	:width: 25%
	:align: center


**Update of Domain bounds** : Once the bounds for the regular bodies and the *local subdomain* has been dispatched, information on the other subdomain bounds are obtained via the function :yref:`yade.mpy.updateDomainBounds`. In this collective communication, each subdomain broadcasts 
its :yref:`Aabb.min` and :yref:`Aabb.max` to other subdomains. Figure `fig-subdomain-bounds`_  shows a schematic in which each subdomain has received the :yref:`Aabb.min` and :yref:`Aabb.max` of the other subdomains. 

.. _fig-subdomain-bounds:
.. figure:: fig/mpy_schema2.*
    :width: 40%
    :align: center
    
**Parallel Collision detection** : 

- Once the  :yref:`Aabb.min` and :yref:`Aabb.max` of the other subdomains are obtained, the collision detection algorithm is used to determine the bodies that have intersections with the remote subdomains. The ids of the identified bodies are 
  then used to build the :yref:`Subdomain.intersections` list. 

 .. _fig-schema-localIntersections:
 .. figure:: fig/mpy_schema3.*
    :width: 40%
    :align: center

- Next step involves in obtaining the ids of the remote bodies intersecting with the current subdomain (:yref:`Subdomain.mirrorIntersections`). Each subdomain sends its list of local body intersections to the respective remote subdomains and also receives the list of intersecting ids from the other subdomains. 
  If the remote bodies do not exist within the current subdomain's :yref:`BodyContainer`, the subdomain then *requests* these remote bodies from the respective subdomain.  A schematic of this operation is shown in figure `fig-schema-mirrorIntersections`_, 
  in which subdomain=1 receives three bodies from subdomain=2, and 1 body from subdomain=0. subdomain=2 receives three bodies from subdomain=1. subdomain=0 only sends its bodies and does *not* receieve from the worker 
  subdomains. This operation sets the stage for communication of the body states to/from the other subdomains. 

 .. _fig-mirrorIntersections:
 .. figure:: fig/mpy_sendBodies.*
    :width: 40%
    :align: center


**Update states** :  

Once the subdomains and the associated intersecting bodies, and remote bodies are identified, :yref:`State` of these bodies are sent and received every timestep, by peer-to-peer communications between the interacting subdomains. In the case of an interaction with the master subdomain (subdomain=0), only the total force and torque exerted on master's bodies by a given subdomain are sent. Figure `fig-sendRecvStates`_ shows a schematic in which the states of the remote bodies between subdomain=1 and subdomain=2 are communicated. Subdomain=0 receives forces and torques from subdomain=1 and subdomain=2. 

.. _fig-sendRecvStates:
.. figure:: fig/mpy_schema4.*
    :width: 40%
    :align: center
        

        
Execution
_________

This section presents methods to execute yade with MPI multiprocessing. In principle the number of processes $np$ can be larger than the number of available cores without problem (this is called oversubscribing, it may also fail depending on OS and MPI implementation). There is no performance gain to expect from oversubscribing, and in production it should be avoided. However it can be useful for experiments (e.g. for testing the examples in this page on a single-core machine).


Interactive mode
----------------
The interactive mode aims primarily at inspecting the simulation after some MPI execution, for debugging for instance. However, functions shown here (especially sendCommand()) may also be usefull to achieve advanced tasks such as controlling transitions between phases of a simulation, collecting and processing results.
The first two flavors may not be used very often in practice, however understanding them is a good way to understand what happens behind the scene.

**Explicit initialization from python prompt**

A pool of yade instances can be spawned with mpy.initialize() as illustrated hereafter. Mind that the next sequences of commands are supposed to be typed directly in the python prompt after starting yade normally, it will not give exactly the same result if it is pasted into a script executed by Yade (see the next section on automatic initialization).

.. initialize the context for next "ipython" sections
.. ipython::
	:suppress:

	Yade [0]: O.reset()

	Yade [1]: from yade.utils import *


.. ipython::
	:verbatim:

	Yade [2]: wallId=O.bodies.append(box(center=(0,0,0),extents=(2,0,1),fixed=True))

	Yade [3]: for x in range(-1,2):
	   ...:    O.bodies.append(sphere((x,0.5,0),0.5))
	   ...:

	Yade [5]: from yade import mpy as mp

	Yade [6]: mp.initialize(3)
	Master: I will spawn  2  workers
	->  [6]: (0, 3)


.. ipython::

	@doctest
	Yade [1]: 1+1
	->  [1]: 4


CODE

After mp.initialize(np) the parent instance of yade takes the role of master process (rank=0). It is the only one executing the commands typed directly in the prompt.
The other instances (rank=1 to rank=np-1) are idle and they wait for commands sent from master.

CODE

Sending commands to the other instances can be done with mpy.sendCommand(), which by default returns the result or the list of results.

CODE (check that scene pointers are different)
CODE (len(bodies) = 1,0,0,0,...)

Sending commands makes it possible to manage all types of message passing using calls to mpi4py. Every picklable python object (namely, nearly all Yade objects) can be transmitted this way:

CODE (send body)
CODE (len(bodies) = 1,0,0,0,...)


**Explicit initialization from python script**

Though usefull, the function sendCommand() is not enough to efficiently manipulate the yade instances in all cases. Even basic features of the python language are missing, e.g. function definitions and loops are a problem - in fact every code fragment which can't fit on a single line is. That is a reason why the mpy module provides a mechanism to initialize from a script.

Whenever Yade is started with a script as argument the script name will be remembered, and if initialize() is executed (in the script itself or interactively in the prompt) all Yade instances will be initialized with that same script. It makes distributing function definitions and simulation parameters trivial (and even distributing scene constructions as seen later). This behaviour is very close to what happens very classicaly in the passive mode, i.e. all processes execute the same program.

If the previous commands are pasted into a script used to start Yade, there is a small surprise, now all instances insert the body.

CODE

That's because all instances executed the script in the initialize() phase. Though logical, this result is not what we want usually if we try to split a simulation into pieces. The solution (typical of all mpi programs) is to use rank of the process in conditionals:

CODE

**Automatic initialization**

Effectively running DEM in parallel on the basis of just the above commands is probably accessible to good hackers but it would be tedious and computationaly innefficient. mpy provides the function mpirun which automatizes most of the steps required for the consistent time integration of a distributed scene, as described in :ref:`introduction <sect_implementation_example2D>`. This includes, mainly, splitting the scene in subdomains based on rank assigned to bodies and handling collisions between the subdomains as time integration proceeds. 

If needed the first execution of mpirun will call the function initialize(), which can therefore be omitted on user's side in most cases.

Here is a concrete example where a floor is assigned to master and multiple groups of spheres are assigned to subdomains:


[CODE] test3D
[COMMENTS] merge/not, erase/master/not, w_interaction/not...


If withMerge=True the bodies in master are updated to reflect in the master scene the evolution of their distributed counterparts. This is done once after finishing the required number of iterations in mpirun. This *merge* operation can include updating interactions.
Merging is an expensive task which requires the communication of large messages and, therefore, it should be done purposely and at a reasonable frequency. It can even be the main bottleneck for massively parallel scenes. Nevertheless it can be usefull for debugging using the 3D view, or for various post-processing tasks. Beyond that it is not required for a proper time integration in general.

**Don't know how to split? Leave it to mpirun**

 mpirun will decide by itself how to distribute the bodies across several subdomains if XXX=True. In such case the difference between the sequential script and its mpi version is limited to importing mpy and calling mpirun after turning the :code:`DOMAIN_DECOMPOSITION` flag ON.  
 
 The automatic splitting of bodies to subdomains is based on the Orthogonal Recursive Bisection Algortithm of Berger [Berger1987]_, and [Fleissner2007]_. The partitioning is based on bisecting the space into several levels, with the longest axis in each level chosen as 
 the bisection axis. The number of levels is determined as :math:`int(log_{2}(N_{w}))` with :math:`N_{w}` being the number of worker subdomains. A schematic of this decomposition is shown in `fig-bisectionAlgo`_, totally there are 4 subdomains. At the initial stage (level = 0),  we assume 
 that subdomain=1 contains the information of the body positions, the longest axis is first determined, this forms the splitting axis/plane. The median of this sorted list is determined, and bodies with positions less than the median is coloured with the current subdomain, (SD=1) and the other half is coloured with 
 SD = 2, the subdomain colouring is determined using the following rule::
      
      if (subdomain <  1<<level) : this subdomain gets the bodies with position lower than the median. 
      if ((subdomain >  1<<level) and (subdomain <  1<<(level+1) ) ) : this subdomain gets the bodies with position greater than median, from subdomain - (1<<level) 
      
     
 This process is continued until the number of levels are reached.
   
 .. _fig-bisectionAlgo:
 .. figure:: fig/mpy_recursuveBisection.*
    :width: 40%
    :align: center

 Figure `fig-domainDecompose`_ shows the resulting partitioning obtained using the ORB algorithm : (a) for 4 subdomains, (b) for 8 subdomains. Odd number of worker subdomains are also supported with the present implementation.
 
 .. _fig-domainDecompose:
 .. figure:: fig/mpy_ddcmp.*
    :width: 40%
    :align: center

 The present implementation can be found in :yref:`yade.bisectionDecomposition`, and a parallel version can be found `here. <https://github.com/bchareyre/yade-mpi/blob/593a4d6abf7e488ab1ac633a1e6725ac301b2a14/py/tree_decomp.py>`_
    


Passive mode
------------

Running in passive mode is straightforward, one just needs to set the number of timesteps as an argument for the :yref:`yade.mpy.mpirun` function. If a scene merge is required, the *withMerge* argument of :yref:`yade.mpy.mpirun` has to be set to true. 
The simulation (:ysrc:`examples/mpi/vtkRecorderExample.py`) is executed with the following command::
  
  mpiexec -np NUMSUBD+1 yade-mpi vtkRecorderExample.py 

where *NUMSUBD* corresponds to the required number of worker subdomains.    


Centralized scene construction
------------------------------
In the centralized method of scene construction, the master subdoamin generates/creates all the bodies (including subdomain), assigns subdomains to tbe bodies and performs the necessary modifications to the engines and collider settings. This scene is  then broadcasted to the workers. 
The workers receives the scene, identifies its respective bodies via :yref:`Body.subdomain` attribute, as worker :code:`rank==b.subdomain` . For large numer of bodies and processes, the centralized scene construction and distribution would take siginificant time for initialization. 


Distributed scene construction
------------------------------

As mentioned in the previous section, the main draw back of the method of centralized scene construction and broadcast leads to long initialization times, an alternative method would be to use the distributed scene construction. 
In this mode of scene construction ther workers first initialize their :yref:`BodyContainer` with the global total number of bodies in the simulation. Each subdomain then creates and inserts bodies at specific location of the initialized but empty :yref:`BodyContainer` using 
:yref:`BodyContainer.insertAtId` function. The user is in charge of setting up the subdomains and partitioning the bodies, an example showing the use of distributed insertion can be found here :ysrc:`examples/mpi/parallelBodyInsert3D.py`. 


Problems to expect
------------------

.. _sect_mpi_reduction

Reduction (partial sums)
------------------------


Control variables
_________________

 - VERBOSE_OUTPUT : Details on each *operation/step* (such as :yref:`yade.mpy.splitScene`, :yref:`yade.mpy.parallelCollide` etc) is printed on the console, useful for debugging purposes
 - ACCUMULATE_FORCES : Control force summation on bodies owned by the master. 
 - ERASE_REMOTE_MASTER : Erase remote bodies in the master subdomain or keep them as unbounded ? Useful for fast merge.
 - OPTIMIZE_COM, USE_CPP_MPI : Use optimized communication functions and MPI functions from :yref:`Subdomain` class 
 - YADE_TIMING : Report timing statistics, prints time spent in communications, collision detection and other operations. 
 - DISTRIBUTED_INSERT : Bodies are created and inserted by each subdomain. 
 - DOMAIN_DECOMPOSITION : If true, the bisection decomposition algorithm is used to assign bodies to the workers/subdomains. 
 - REALLOCATE_FREQUENCY : if > 0, bodies are migrated between subdomains for efficient load balancing. 
 - FLUID_COUPLING : Flag for coupling with OpenFOAM. 
 


Various remarks
_______________
- sendCommand() has a hardcoded latency of 0.001s to not keep all cores 100\% busy waiting for a command (with possibly little left to OS). If sendCommand() is used at high frequency in complex algorithms it might be beneficial to decrease that sleep time.
