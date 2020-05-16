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

**subdomain**: a (sub)set of bodies attached to one MPI process after domain decomposition - with or without spatial coherence. The corresponding class in Yade is :yref:`Subdomain`, a `Shape` instance with helper functions for MPI communications. In some sense `Subdomain` is to subscribed bodies what :yref:`Clump` (another `Shape`) is to clump members.

**rank**: subdomain index from 0 to *N*-1  (with *N* the number of mpi processes) to identify subdomains. The rank of the subdomain a body belongs to can be retrieved as :yref:`Body.subdomain`. Each subdomain corresponds to an instance of yade and a specific scene during parallel execution. The rank of the scene is given by :yref:`Scene.subdomain`.

**master**: refers to subdomain with *rank* =0. This subdomain does not behave like others. In general master will handle boundary conditions and it will control transitions and termination of the whole simulation. Unlike standard subdomains it may not contain a large number of raw bodies (i.e. not beyond objects bounding the scene such as walls or boxes). In interactive execution master is the process responding to the python prompt.

**splitting and merging**: cutting a master :yref:`Scene` into a set of smaller, distributed, scenes is called "splitting". The split is undone by a 'merge', by which all bodies and (optionally) all interactions are sent back to the master thread. Splitting, running, then merging, should leave the scene just as if no MPI had been used at all (i.e. as if the same number of iterations had been executed in single-thread). Therefore normal O.run() after that should work as usual.

**intersections**: subsets of bodies in a subdomain intersected by the bounding box of other subdomains (see `fig-subdomains`_). *intersection(i,j)* refers to the bodies owned by current (*i*) subdomain and intersecting subdomain *j* (retrieved as :yref:`O._sceneObj.subD.intersections[j]<Subdomain.intersections>`); *mirrorIntersection(i,j)* refers to bodies owned by *j* and intersecting current domain (retrieved as :yref:`O._sceneObj.subD.mirrorIntersections[j]<Subdomain.mirrorIntersections>`). The bodies are listed by :yref:`Body.id`. By definition *intersection(i,j)=mirrorIntersection(j,i)*.

The intersections and mirror intersections are updated automatically as part of parallel collision detection. They define which body states need to be communicated. The bodies in intersections need to be *sent* to other subdomains (in pratice only updated position and velocity are sent at every iteration), the bodies in mirrorIntersections need to be received from other subdomains.


.. _fig-subdomains:
.. figure:: fig/subdomains.*
	:width: 12cm
	:align: center

Two overlapping subdomains and their intersections. In this situation we have *SubD1.intersections[SubD2.subdomain]=[id4,id5]* and *SubD1.mirrorIntersections[SubD2.subdomain]=[id1]*, with *SubD1* and *SubD2* instances of :yref:`Subdomain`.


.. _sect_mpi_implementation:

Walkthrough
-----------

For demonstrating the main internal steps in the implemented parallel algorithm let us consider the example script :ysrc:`examples/mpi/testMPI_2D.py`. Executing this script (interactive or passive mode) with three MPI processes generates the scene as shown in `fig-scene`_. It then executes :yref:`mpirun<yade.mpy.mpirun>`, which triggers the steps described hereafter.

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

- Next step involves obtaining the ids of the remote bodies intersecting with the current subdomain (:yref:`Subdomain.mirrorIntersections`). Each subdomain sends its list of local body intersections to the respective remote subdomains and also receives the list of intersecting ids from the other subdomains. 
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

This section presents methods to execute yade with MPI support. The mpy modules tries to retain from yade design an important feature: interactive access to the objects of scene (or of multiple scenes in this case), as explained below. Interactive execution does not use the `mpiexec` command of OpenMPI, a pool of workers is spawned by the mpy module after yade startup. In production one may use passive jobs instead, and in that case `mpiexec` will preceed the call to yade.

.. note:: Most examples in this page use 4 mpi processes. It is not a problem, in principle, to run the examples even if the number of available cores is less than 4 (this is called oversubscribing, it may also fail depending on OS and MPI implementation). There is no performance gain to expect from oversubscribing, and in production it should be avoided, but it is useful for experiments (e.g. for testing the examples in this page on a single-core machine).


Interactive mode
----------------
The interactive mode aims primarily at inspecting the simulation after some MPI execution for debugging. However, functions shown here (especially `sendCommand`) may also be usefull to achieve advanced tasks such as controlling transitions between phases of a simulation, collecting and processing results.
The first two flavors may not be used very often in practice, however understanding them is a good way to understand what happens behind the scene.

**Explicit initialization from python prompt**

A pool of yade instances can be spawned with mpy.initialize() as illustrated hereafter. Mind that the next sequences of commands are supposed to be typed directly in the python prompt after starting yade, it will not give exactly the same result if it is pasted into a script executed by Yade (see the next section on automatic initialization).

.. initialize the context for next "ipython" sections
.. ipython::
	:suppress:

	Yade [0]: O.reset()

	Yade [1]: from yade.utils import *
	
.. Note: most of the blocks in next ipython directives executes correctly, however it gives slightly different output through sphinx pipeline, unfortunately. Comment out the verbatim to make sphinx effectively execute these mpi sections. 
	
.. ipython::
	:verbatim:

	Yade [2]: wallId=O.bodies.append(box(center=(0,0,0),extents=(2,0,1),fixed=True))

	Yade [3]: for x in range(-1,2):
	   ...:    O.bodies.append(sphere((x,0.5,0),0.5))
	   ...:

	Yade [5]: from yade import mpy as mp

	@doctest
	Yade [6]: mp.initialize(4)
	Master: I will spawn  3  workers
	->  [6]: (0, 4)


.. .. ipython::
.. 
.. 	@doctest
.. 	Yade [1]: 1+1
.. 	->  [1]: 4


After mp.initialize(np) the parent instance of yade takes the role of master process (rank=0). It is the only one executing the commands typed directly in the prompt.
The other instances (rank=1 to rank=np-1) are idle and they wait for commands sent from master. Sending commands to the other instances can be done with `mpy.sendCommand()`, which by default returns the result or the list of results. We use that command below to verify that the spawned workers point to different (still empty) scenes:

.. ipython::
	:verbatim:
	
	Yade [8]: len(O.bodies)
	 ->  [8]: 4

	Yade [9]: mp.sendCommand(executors="all",command="str(O)") # check scene pointers
	->  [9]: ['<yade.wrapper.Omega object at 0x7f6db7012300>', '<yade.wrapper.Omega object at 0x7f94c79ec300>', '<yade.wrapper.Omega object at 0x7f5519742300>', '<yade.wrapper.Omega object at 0x7f264dd80300>']

	Yade [10]: mp.sendCommand(executors="all",command="len(O.bodies)",wait=True) #check content
	->  [10]: [4, 0, 0, 0]

Sending commands makes it possible to manage all types of message passing using calls to the underlying mpi4py (see mpi4py documentation for more functionalities).

.. ipython::
	:verbatim:
	
	Yade [3]: mp.sendCommand(executors=1,command="message=comm.recv(source=0); print('received',message)")

	Yade [4]: mp.comm.send("hello",dest=1)
	received hello

Every picklable python object (namely, nearly all Yade objects) can be transmitted this way. Remark hereafter the use of :yref:`mpy.mprint <yade.mpy.mprint>` (identifies the worker by number and by font colors). Note also that the commands passed via `sendCommand` are executed in the context of the mpy module, for this reason `comm`, `mprint`, `rank` and all objects of the module are accessed without the `mp.` prefix.

.. ipython::
	:verbatim:

	Yade [3]: mp.sendCommand(executors=1,command="O.bodies.append(comm.recv(source=0))",wait=False) # leaves the worker idle waiting for an argument to append()

	Yade [4]: b=Body(shape=Sphere(radius=0.7))  # now create body in the context of master

	Yade [5]: mp.comm.send(b,dest=1) # send it to worker 1

	Yade [6]: mp.sendCommand(executors="all",command="mprint('received',[b.shape.radius if hasattr(b.shape,'radius') else None for b in O.bodies])")
	Master: received [None, 0.5, 0.5, 0.5] 
	Worker1: received [0.7] 
	Worker3: received [] 
	Worker2: received [] 
	->  [5]: [None, None, None, None] # printing yields no return value, hence that empty list, "wait=False" argument to sendCommand would suppress it


**Explicit initialization from python script**

Though usefull for advanced operations, the function sendCommand() is limited. Basic features of the python language are missing, e.g. function definitions and loops are a problem - in fact every code fragment which can't fit on a single line is. In practice the mpy module provides a mechanism to initialize from a script, where functions and variables will be declared.

Whenever Yade is started with a script as argument the script name will be remembered, and if mpy.initialize() is called (by the script itself or interactively in the prompt) all Yade instances will be initialized with that same script. It makes distributing function definitions and simulation parameters trivial (and even distributing scene constructions as seen below). This behaviour is what happens classicaly with MPI: all processes execute the same program.

If the first commands above are pasted into a script used to start Yade, there is a small surprise: all instances insert the same bodies as master (with interactive execution only master was inserting). Here is the script:

.. code-block:: python

	# script 'test1.py'
	wallId=O.bodies.append(box(center=(0,0,0),extents=(2,0,1),fixed=True))
	for x in range(-1,2):
		O.bodies.append(sphere((x,0.5,0),0.5))
	from yade import mpy as mp
	mp.initialize(4)
	print( mp.sendCommand(executors="all",command="str(O)",wait=True) )
	print( mp.sendCommand(executors="all",command="len(O.bodies)",wait=True) )

and the output reads:

.. code-block:: none
	
	yade test1.py 
	...
	Running script test1.py
	Master: will spawn  3  workers 
	None
	None
	None
	None
	None
	None
	['<yade.wrapper.Omega object at 0x7feb979403a0>', '<yade.wrapper.Omega object at 0x7f5b61ae9440>', '<yade.wrapper.Omega object at 0x7fdd466b8440>', '<yade.wrapper.Omega object at 0x7f8dc7b73440>']
	[4, 4, 4, 4]

That's because all instances executed the script in the initialize() phase. "None" is printed 2x3 times because the script contains `print( mp.sendCommand(...))` twice, the workers try to execute that too, but for them `sendCommand` returns by default, hence the None.


Though logical, this result is not what we want usually if we try to split a simulation into pieces. The solution (typical of all mpi programs) is to use the `rank` of the process in conditionals. In order to produce the same result as before, for instance, the script can be modified as follows.

.. code-block:: python

	# script 'test2.py'
	from yade import mpy as mp
	mp.initialize(4)
	if mp.rank==0: # only master
		wallId=O.bodies.append(box(center=(0,0,0),extents=(2,0,1),fixed=True))
		for x in range(-1,2):
		O.bodies.append(sphere((x,0.5,0),0.5))

		print( mp.sendCommand(executors="all",command="str(O)",wait=True) )
		print( mp.sendCommand(executors="all",command="len(O.bodies)",wait=True) )
		print( mp.sendCommand(executors="all",command="str(O)",wait=True) )
		
.. code-block:: none

	Running script test2.py
	Master: will spawn  3  workers 
	['<yade.wrapper.Omega object at 0x7f21a8c8d3a0>', '<yade.wrapper.Omega object at 0x7f3142e43440>', '<yade.wrapper.Omega object at 0x7fb699b1a440>', '<yade.wrapper.Omega object at 0x7f1e4231e440>']
	[4, 0, 0, 0]


We could also use `rank` to assign bodies from different regions of space to different workers, as found in example :ysrc:`examples/mpi/helloMPI.py`, with rank-dependent positions:

.. code-block:: python

	# rank is accessed without "mp." prefix as it is interpreted in mpy module's scope
	mp.sendCommand(executors=[1,2],command= "ids=O.bodies.append([sphere((xx,1.5+rank,0),0.5) for xx in range(-1,2)])")
	


	
**Automatic initialization**

Effectively running a distibuted DEM simulation on the basis of just the above commands would still be tedious and computationaly innefficient. The mpy modules thus provides the function :yref:`mpy.mpirun <yade.mpy.mpirun>` to automatize most of the steps, as described in :ref:`introduction <sect_implementation_example2D>`. This includes, mainly, splitting the scene in subdomains based on rank assigned to bodies and handling collisions between the subdomains as time integration proceeds. 

If needed the first execution of mpirun will call the function initialize(), which can therefore be omitted on user's side in most cases.

Here is a concrete example where a floor is assigned to master and multiple groups of spheres are assigned to subdomains::

	NSTEPS=5000 #turn it >0 to see time iterations, else only initilization 
	numThreads = 4 # number of threads to be spawned, (in interactive mode).

	import os
	from yade import mpy as mp

	#materials 
	young = 5e6
	compFricDegree = 0.0
	O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle = radians(compFricDegree), density= 2600, label='sphereMat'))
	O.materials.append(FrictMat(young=young*100, poisson = 0.5, frictionAngle = compFricDegree, density =2600, label='wallMat'))


	#add spheres
	
	mn,mx=Vector3(0,0,0),Vector3(90,180,90)
	pred = pack.inAlignedBox(mn,mx)
	O.bodies.append(pack.regularHexa(pred,radius=2.80,gap=0, material='sphereMat'))

	#walls (floor)
	
	wallIds=aabbWalls([Vector3(-360,-1,-360),Vector3(360,360,360)],thickness=10.0, material='wallMat')
	O.bodies.append(wallIds)

	#engines 
	O.engines=[
		ForceResetter(),
		InsertionSortCollider([
			Bo1_Sphere_Aabb(),
			Bo1_Box_Aabb()], label = 'collider'), # always add labels. 
		InteractionLoop(
			[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
			[Ip2_FrictMat_FrictMat_FrictPhys()],
			[Law2_ScGeom_FrictPhys_CundallStrack()], 
			label="interactionLoop"
		),
		GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.3,  timeStepUpdateInterval=100, parallelMode=True, label = 'timeStepper'),
		NewtonIntegrator(damping=0.1,gravity = (0, -0.1, 0), label='newton'), 
		VTKRecorder(fileName='spheres/3d-vtk-', recorders=['spheres', 'intr', 'boxes'], parallelMode=True,iterPeriod=500), #use .pvtu to open spheres, .pvtp for ints, and .vtu for boxes.
	]

	#set a custom verletDist for efficiency. 
	collider.verletDist = 1.5

	#########  RUN  ##########
	# customize mpy
	mp.ERASE_REMOTE_MASTER = True   #keep remote bodies in master? 
	mp.DOMAIN_DECOMPOSITION= True	#automatic splitting/domain decomposition
	#mp.mpirun(NSTEPS)		#passive mode run 
	mp.MERGE_W_INTERACTIONS = False
	mp.mpirun(NSTEPS,numThreads,withMerge=True) # interactive run, numThreads is the number of workers to be initialized, see below for withMerge explanation.
	mp.mergeScene()  #merge scene after run. 
	if mp.rank == 0: O.save('mergedScene.yade')

	#demonstrate getting stuff from workers, here we get kinetic energy from worker subdomains, notice that the master (mp.rank = 0), uses the sendCommand to tell workers to compute kineticEnergy. 
	if mp.rank==0:
		print("kinetic energy from workers: "+str(mp.sendCommand([1,2],"kineticEnergy()",True)))
		

The script is then executed as follows::
  
  yade script.py 

For running further timesteps, the mp.mpirun command has to be executed in yade prompt:
  
.. ipython::
	:verbatim:
.. 	:suppress:
	
	Yade [0]: mp.mpirun(100,4,withMerge=False) #run for 100 steps and no scene merge. 
	
	Yade [1]: mp.sendCommand([1,2],"kineticEnergy()",True) # get kineticEnergy from workers 1 and 2. 
	
	Yade [2]: mp.mpirun(1,4,withMerge=True) #run for 1 step and merge scene into master. 
	


**Don't know how to split? Leave it to mpirun**

Splitting an initial scene into subdomains and updating the subdomains after particle motion are critical issues in terms of efficiency.
The mpi implementation has algorithms for both.

.. note:: The mpy module has no requirement in terms of how the subdomains are defined. Even assigning the bodies randomly to subdomains would work. It would only be suboptimal as the number of interactions between subdomains would increase compared to a proper partition of space. Using the helper functions described here is not necessary.


**Initial split**
 mpirun will decide by itself how to distribute the bodies across several subdomains if *DOMAIN_DECOMPOSITION* =True. In such case the difference between the sequential script and its mpi version is limited to importing mpy and calling mpirun after turning the *DOMAIN_DECOMPOSITION* flag.  
 
 The automatic splitting of bodies to subdomains is based on the Orthogonal Recursive Bisection Algortithm of Berger [Berger1987]_, and [Fleissner2007]_. The partitioning is based on bisecting the space at several *levels*, with the longest axis in each level chosen as 
 the bisection axis. The number of levels is determined as :math:`int(log_{2}(N_{w}))` with :math:`N_{w}` being the number of worker subdomains. A schematic of this decomposition is shown in `fig-bisectionAlgo`_, with 4 worker subdomains. At the initial stage (level = 0),  we assume 
 that subdomain=1 contains the information of the body positions (and bodies), the longest axis is first determined, this forms the bisectioning axis/plane. The list containing the body positions is sorted along the bisection axis, and the median of this sorted list is determined. The bodies with positions (bisection coordinate) less than the median is coloured with the current subdomain, (SD=1) and the other half is coloured with 
 SD = 2, the subdomain colouring at each level is determined using the following rule::
      
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

 The present implementation can be found in :ysrc:`py/bisectionDecomposition.py`, and a parallel version can be found `here. <https://github.com/bchareyre/yade-mpi/blob/593a4d6abf7e488ab1ac633a1e6725ac301b2a14/py/tree_decomp.py>`_

 
**Updating the decomposition (dynamic load balancing)**

As the bodies move each subdomain may experience overall distorsion and diffusion of bodies to/from other subdomains. We want to keep the subdomains as compact as possible to minimize communications, instead. An algorithm does that dynamically if :yref:`mpy.REALLOCATE_FREQUENCY <yade.mpy.REALLOCATE_FREQUENCY>`>0. It exploits :yref:`InsertionSortCollider` to reassign bodies efficiently and in synchronicity with collision detection.

The algorithm is *not* centralized in order to preserve scalability. It involves only peer-to-peer communications between the workers which share an intersection.

The criterion for re-allocating bodies involves finding the position a median plane between two subdomains such that after discriminating bodies on the "+" and "-" side of that plane the total number in each subdomain is preserved. It results in the type of split shown in the video hereafter. Even though the median planes seem to rotate rather quickly at some point in this video, there are actually five collision detections between each re-allocation, i.e. thousands of time iterations to effectively rotate the split between two different colors. These progressive rotations are for good since the initial split would have resulted in flat discs, mostly.

.. youtube:: Qb5vPjRPFRw

**Merging**

The possibility of a "merge" emerged in previous example, as an optional argument of `mpirun` or as a standalone function :yref:`mpy.mergeScene <yade.mpy.mergeScene>`. 

If withMerge=True in mpirun then the bodies in master scene are updated to reflect the evolution of their distributed clones. This is done once after finishing the required number of iterations in mpirun. This *merge* operation can include updating interactions. :yref:`mpy.mergeScene <yade.mpy.mergeScene>` does the same within current iteration.
Merging is an expensive task which requires the communication of large messages and, therefore, it should be done purposely and at a reasonable frequency. It can even be the main bottleneck for massively parallel scenes. Nevertheless it can be usefull for debugging using the 3D view, or for various post-processing tasks. 
The *MERGE_W_INTERACTIONS* provides full merge, i.e. the interactions in the worker subdomains and between the subdomains are included, else only the position and states of the bodies are use. Merging with interactions should result in a usual yade scene, ready for further time-stepping in non-mpi mode or (more usefull) for some post-processing. The merge operation is not required for a proper time integration in general.


For MPI cases, the *parallelMode* flag for :yref:`GlobalStiffnessTimeStepper` and :yref:`VTKRecorder` have to be turned on. 

Passive mode
------------

Running in passive mode is straightforward, one just needs to set the number of timesteps as an argument for the :yref:`yade.mpy.mpirun` function. If a scene merge is required, the *withMerge* argument of :yref:`yade.mpy.mpirun` has to be set to true. 
The simulation (:ysrc:`examples/mpi/vtkRecorderExample.py`) is executed with the following command::
  
  mpiexec -np NUMSUBD+1 yade vtkRecorderExample.py 

where *NUMSUBD* corresponds to the required number of worker subdomains.    


Centralized versus distributed scene construction
-------------------------------------------------
In the centralized method of scene construction, the master process creates all the bodies of a scene and assigns subdomains to them. As part of mpy initialization some engines will be modified or inserted, then the scene is broadcasted to the workers.
Each worker receives the entire scene, identifies its assigned bodies via :yref:`Body.subdomain` (if worker's :code:`rank==b.subdomain` the bodies are retained) and erase the others. Such a scene construction was used in the previous example and it is by far the simplest. It makes no real difference with building a scene for non-MPI execution besides calling `mp.mpirun` instead or jusr `O.run`.

For large number of bodies and processes, though, the centralized scene construction and distribution can take a significant time. It can also be memory bound since the memory usage is quadratic: suppose N bodies per thread on a 32-core node, centralized construction implies that 32 copies of the entire scene exist simultaneously in memory at some point in time (during the split), i.e. $32^2 N$ bodies on one single node. For massively parallel applications distributed construction should be prefered.

In distributed mode each worker instantiates its own bodies and insert them in the local :yref:`BodyContainer`. Attention need to be paid to properly assign bodies ids since no index should be owned by two different workers initially. Insertion of bodies  in :yref:`BodyContainer` with imposed ids is done with 
:yref:`BodyContainer.insertAtId`. The distributed mode is activated by setting the :code:`DISTRIBUTED_INSERT` flag ON, the user is in charge of setting up the subdomains and partitioning the bodies, an example showing the use of distributed insertion can be found in :ysrc:`examples/mpi/parallelBodyInsert3D.py`. 

The bissection algorithm can be used for defining the initial split, in the distributed case too, since it takes a points dataset as input. Provided that all workers work with the same dataset (e.g. the same sequence of a random number generator) they will all reach the same partitioning, and they can instanciate their bodies on this basis. 


Problems to expect
------------------

.. _sect_mpi_reduction

Reduction (partial sums)
------------------------

Quantities such as kinetic energy cannot be obtained for the entire scene just by summing the return value of `kineticEnergy()` from each subdomain.
This is because each subdmomain may contain also images of bodies from intersecting subdomains and they may add their velocity, mass, or whatever is summed, to what is returned by each worker.
Although some most-used functions of yade may progressively get mpi support to filter out bodies from remote domains, it is not standard yet and therefore partial sums may need to be implemented on a case-by-case basis.




Control variables
_________________

 - VERBOSE_OUTPUT : Details on each *operation/step* (such as :yref:`yade.mpy.splitScene`, :yref:`yade.mpy.parallelCollide` etc) is printed on the console, useful for debugging purposes
 - ACCUMULATE_FORCES : Control force summation on bodies owned by the master. 
 - ERASE_REMOTE_MASTER : Erase remote bodies in the master subdomain or keep them as unbounded ? Useful for fast merge.
 - OPTIMIZE_COM, USE_CPP_MPI : Use optimized communication functions and MPI functions from :yref:`Subdomain` class 
 - YADE_TIMING : Report timing statistics, prints time spent in communications, collision detection and other operations. 
 - DISTRIBUTED_INSERT : Bodies are created and inserted by each subdomain, used for distributed scene construction. 
 - DOMAIN_DECOMPOSITION : If true, the bisection decomposition algorithm is used to assign bodies to the workers/subdomains. 
 - MINIMAL_INTERSECTIONS : Reduces the size of position/velocity communications (at the end of the colliding phase, we can exclude those bodies with no interactions besides body<->subdomain from intersections). 
 - REALLOCATE_FREQUENCY : if > 0, bodies are migrated between subdomains for efficient load balancing. If =1 realloc. happens each time collider is triggered, else every N collision detection
 - REALLOCATE_MINIMAL : Intersections are minimized before reallocations, hence minimizing the number of reallocated bodies
 - USE_CPP_REALLOC : Use optimized C++ functions to perform body reallocations
 - FLUID_COUPLING : Flag for coupling with OpenFOAM. 
 

Various remarks
_______________
- sendCommand() has a hardcoded latency of 0.001s to not keep all cores 100\% busy waiting for a command (with possibly little left to OS). If sendCommand() is used at high frequency in complex algorithms it might be beneficial to decrease that sleep time.
