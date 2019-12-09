

wallId=O.bodies.append(box(center=(0,0,0),extents=(2,0,1),fixed=True))
for x in range(-1,2):
	O.bodies.append(sphere((x,0.5,0),0.5))

newton.gravity=(0,-10,0) #else nothing would move
O.dt=0.1*PWaveTimeStep() 
O.dynDt=False #very important, we don't want subdomains to use many different timesteps...

from yade import mpy as mp
mp.initialize(3)

if mp.rank==0:
	print( mp.sendCommand(executors="all",command="len(O.bodies)",wait=True) )

	#note that 'rank' is used instead of mp.rank, the scope of the mpy module is accessed here
	mp.sendCommand(executors=[1,2],command= "ids=O.bodies.append([sphere((xx,1.5+rank,0),0.5) for xx in range(-1,2)])",wait=True)

	print( mp.sendCommand(executors="all",command="len(O.bodies)",wait=True) )

	# b.subdomain is defined for all bodies; 'for ...' loops are awkward in a one-liner hence 'map' is used instead
	mp.sendCommand(executors=[1,2],command= "list(map(lambda b: setattr(b,'subdomain',rank),O.bodies))", wait=True)

	print("Assigned bodies:", mp.sendCommand([1,2], "len([b for b in O.bodies if b.subdomain==rank])", True) )
	
#mp.DISTRIBUTED_INSERT=True
mp.MERGE_W_INTERACTIONS=True
mp.ERASE_REMOTE_MASTER=False

# we can exploit sendCommand to send data between yade instances directly with mpi4py functions (see mpi4py documentation)
#Yade [8]: mp.sendCommand(executors=1,command="message=comm.recv(source=0); print('received',message)",wait=False)

#Yade [10]: mp.comm.send("hello",dest=1)
#received hello

#Yade [11]: mp.sendCommand(executors=1,command="message=comm.recv(source=0); mprint('received: ',message)",wait=False)

#Yade [12]: mp.comm.send("hello",dest=1)
#Worker1: received hello


#inspect 3D view



##########  RUN  ##########
#def collectTiming():
	#created = os.path.isfile("collect.dat")
	#f=open('collect.dat','a')
	#if not created: f.write("numThreads mpi omp Nspheres N M runtime \n")
	#from yade import timing
	#f.write(str(numThreads)+" "+str(os.getenv('OMPI_COMM_WORLD_SIZE'))+" "+os.getenv('OMP_NUM_THREADS')+" "+str(N*M*(numThreads-1))+" "+str(N)+" "+str(M)+" "+str(timing.runtime())+"\n")
	#f.close()


#if rank is None: #######  Single-core  ######
	#O.timingEnabled=True
	#O.run(NSTEPS,True)
	##print "num bodies:",len(O.bodies)
	#from yade import timing
	#timing.stats()
	#collectTiming()
	#print("num. bodies:",len([b for b in O.bodies]),len(O.bodies))
	#print("Total force on floor=",O.forces.f(WALL_ID)[1])
#else: #######  MPI  ######
	## customize
	#mp.ACCUMULATE_FORCES=True #trigger force summation on master's body (here WALL_ID)
	#mp.VERBOSE_OUTPUT=False
	#mp.ERASE_REMOTE=False #erase bodies not interacting wit a given subdomain?
	#mp.OPTIMIZE_COM=True #L1-optimization: pass a list of double instead of a list of states
	#mp.USE_CPP_MPI=True and mp.OPTIMIZE_COM #L2-optimization: workaround python by passing a vector<double> at the c++ level
	#mp.MERGE_SPLIT=mergeSplit
	#mp.COPY_MIRROR_BODIES_WHEN_COLLIDE = bodyCopy and not mergeSplit

	#mp.mpirun(NSTEPS)
	#print ("num. bodies:",len([b for b in O.bodies]),len(O.bodies))
	#if rank==0:
		#mp.mprint( "Total force on floor="+str(O.forces.f(WALL_ID)[1]))
		#collectTiming()
	#else: mp.mprint( "Partial force on floor="+str(O.forces.f(WALL_ID)[1]))
	##mp.mergeScene()
	#if rank==0: O.save('mergedScene.yade')
	##mp.MPI.Finalize()
##exit()

