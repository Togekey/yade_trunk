

import os
from yade import mpy as mp
rank,numThreads = mp.initialize(4)

wallId=O.bodies.append(box(center=(numThreads*N*0.5,-0.5,0),extents=(2*numThreads*N,0,2),fixed=True))
for x in range(3)
	O.bodies.append(sphere((0.5+x,0.5,0),0.5))

newton.gravity=(0,-10,0) #else nothing would move
O.dt=0.1*PWaveTimeStep() #very important, we don't want subdomains to use many different timesteps...
O.dynDt=False


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

