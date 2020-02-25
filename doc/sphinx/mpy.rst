.. _mpy:

MPI parallelization
===================

The module mpy :yref:`yade.mpy` implements parallelization by domain decomposition (distributed memory) using the Message Passing Interface (MPI) implemented by OpenMPI. It aims at exploiting large numbers of cores, where shared memory techniques (OpenMP) are helpess.
The shared memory and the distributed memory approaches are compatible, i.e. it is possible to run hybrid jobs using both, as shown below (and it may well be the optimal solution in many cases).

Most calls to OpenMPI library are done in Python using mpi4py. For the sake of efficiency some critical communications are triggered via python wrappers of C++ functions, wherein messages are produced, sent/received, and processed.

Disclaimer: even though the mpy module provides the function :yref:`yade.mpy.mpirun`, which may seem as a simple replacement for O.run(), setting up a simulation with mpy might be deceptively triavial.
It is anticipated that, in general, a simple replacement of "run" by "mpirun" in an arbitrary script will not work. To understand why, and to tackle the problems, basic knowledge of how MPI works will certainly help (specifically mpi4py).



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

Implementation
---------------

For demonstrating the steps in the implemented parallel algorithm let us conider the example script :ysrc:`examples/mpi/testMPI_2D.py`. Executing this script (interactive or passive mode) with three MPI processes generates the following scene as shown in `fig-scene`_. 

.. _fig-scene:
.. figure:: fig/schema0.*
	:width: 9cm
	:align: center


In this scene, we have three MPI processes (three subdomains)  and the raw bodies are partitioned among the subdomains/ranks 1 and 2. The master process with subdomain=0 holds the boundary/wall type body. Bodies can be manually assigned or automatically assigned via a domain decomposition algorithm. Details 
on the dommain decomposition algorithm is presented in the later section of this documnent. 

In the function :yref:`yade.mpy.splitScene()`, the collider settings are changed, in particular the verlet sweep distance of the :yref:`InsertionSortCollider.sweepLength` is extended. The verletDist is an important factor controlling the efficiency of the simulations. The reason for this will become evident in the later steps. 

**Bounds dispatching** : In the next step, the :yref:`Body.bound` is dispatched with the :yref:`Aabb` extended as shown in figure `fig-regularbounds`_ (in dotted lines). Note that the :yref:`Subdomain` :yref:`Aabb` is obtained from taking the min and max of the owned bodies, see second half of 
figure `fig-regularbounds`_, with solid coloured lines for the subdomain :yref:`Aabb`. At this time, the min and max of other subdomains are unknown. 

.. _fig-regularbounds:
.. figure:: fig/schema1.*
	:width: 9cm
	:align: center

**Update of Domain bounds** : Once the bounds for the regular bodies and the *local subdomain* has been dispatched, information on the other subdomain bounds are obtained via the function :yref:`yade.mpy.updateDomainBounds`. In this collective communication, each subdomain broadcasts 
its :yref:`Aabb.min` and :yref:`Aabb.max` to other subdomains. Figure `fig-subdomain-bounds`_  shows a schematic in which each subdomain has received the :yref:`Aabb.min` and :yref:`Aabb.max` of the other subdomains. 

.. _fig-subdomain-bounds:
.. figure:: fig/schema2.*
    :width: 16cm
    :align: center
    
**Parallel Collision detection** : 

- Once the  :yref:`Aabb.min` and :yref:`Aabb.max` of the other subdomains are obtained, the collision detection algorithm is used to determine the bodies that have intersections with the remote subdomains. The ids of the identified bodies are 
  then used to build the :yref:`Subdomain.intersections` list. 

 .. _fig-schema-localIntersections:
 .. figure:: fig/schema3.*
    :width: 16cm
    :align: center

- Next step involves in obtaining the ids of the remote bodies intersecting with the current subdomain (:yref:`Subdomain.mirrorIntersections`). Each subdomain sends its list of local body intersections to the respective remote subdomains and also receives the list  of intersecting ids from theother subdomains. 
  If the remote bodies do not exist within the current subdomain's :yref:`BodyContainer`, the subdomain then *requests* these remote bodies from the respective subdomain.  A schematic of this operation is shown in figure `fig-schema-mirrorIntersections`_, 
  in which subdomain=1 receives three bodies from subdomain=2, and 1 body from subdomain=0. subdomain=2 receives three bodies from subdomain=1. subdomain=0 only sends its bodies and does NOT receieve from the worker 
  subdomains. This operation sets the stage for communication of the body states to/from the other subdomains. 

 .. _fig-mirrorIntersections:
 .. figure:: fig/sendbodies.*
    :width: 16cm
    :align: center


**Update states** :  

Once the subdomains and the associated intersecting bodies, and remote bodies are identified, :yref:`State` of these bodies are sent and received every timestep. In the case of an interaction with the master subdomain (subdomain=0), the interaction forces and torques
are sent. Figure `fig-sendRecvStates`_ shows a schematic in which the states of the remote bodies between subdomain=1 and subdomain=2 are communicated. Subdomain=0 receives forces and torques from subdomain=1 and subdomain=2. 

.. _fig-sendRecvStates:
.. figure:: fig/schema4.*
    :width: 16cm
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

Effectively running DEM in parallel on the basis of just the above commands is probably accessible to good hackers but it would be tedious and computationaly innefficient. mpy provides the function mpirun which automatizes most of the steps required for the consistent time integration of a distributed scene. This includes, mainly, splitting the scene in subdomains based on indices assigned to bodies and handling collisions between the subdomains as time integration proceeds.

If needed the first execution of mpirun will call the function initialize(), which can therefore be omitted on user's side in most cases.

Here is a concrete example where a floor is assigned to master and multiple groups of spheres are assigned to subdomains:


[CODE] test3D
[COMMENTS] merge/not, erase/master/not, w_interaction/not...


If withMerge=True the bodies in master are updated to reflect in the master scene the evolution of their distributed counterparts. This is done once after finishing the required number of iterations in mpirun. This *merge* operation can include updating interactions.
Merging is an expensive task which requires the communication of large messages and, therefore, it should be done purposely and at a reasonable frequency. It can even be the main bottleneck for massively parallel scenes. Nevertheless it can be usefull for debugging using the 3D view, or for various post-processing tasks. Beyond that it is not required for a proper time integration in general.

**Don't know how to split? Leave it to mpirun**

 mpirun will decide by itself how to distribute the bodies across several subdomains if XXX=True. In such case the difference between the sequential script and its mpi version is limited to importing mpy and calling mpirun after turning that flag on.

 [CODE]
 [BRIEF NOTES ON BISSECTION ALGORITHM - reference?]
 
 .. _fig-bisectionAlgo:
 .. figure:: fig/recursive_bisection.*
    :width: 16cm
    :align: center


 .. _fig-domainDecompose:
 .. figure:: fig/dd1.*
    :width: 16cm
    :align: center



Passive mode
------------





Centralized scene construction
------------------------------

Distributed scene construction
------------------------------

Problems to expect
------------------

Reduction (partial sums)


Control variables
_________________

 - VERBOSE_OUTPUT


Various remarks
_______________
- sendCommand() has a hardcoded latency of 0.001s to not keep all cores 100\% busy waiting for a command (with possibly little left to OS). If sendCommand() is used at high frequency in complex algorithms it might be beneficial to decrease that sleep time.
