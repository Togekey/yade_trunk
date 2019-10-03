
'''
# Possible executions of this script
### Parallel:
# mpiexec -n 4 yade-mpi -n -x testMPIxNxM.py
# mpiexec -n 4 yade-mpi  -n -x testMPIxN.py N M # (n-1) subdomains with NxM spheres each
'''


NSTEPS=1000 #turn it >0 to see time iterations, else only initilization TODO!HACK
#NSTEPS=50 #turn it >0 to see time iterations, else only initilization
N=50; M=50; #(columns, rows) per thread

#################
# Check MPI world
# This is to know if it was run with or without mpiexec (see preamble of this script)
import os
rank = os.getenv('OMPI_COMM_WORLD_RANK')
if rank is not None: #mpiexec was used
	rank=int(rank)
	numThreads=int(os.getenv('OMPI_COMM_WORLD_SIZE'))
else: #non-mpi execution, numThreads will still be used as multiplier for the problem size (2 => multiplier is 1)
	numThreads=2 if len(sys.argv)<4 else (int(sys.argv[3]))
	print("numThreads",numThreads)
	
if len(sys.argv)>1: #we then assume N,M are provided as 1st and 2nd cmd line arguments
	N=int(sys.argv[1]); M=int(sys.argv[2])

############  Build a scene (we use Yade's pre-filled scene)  ############


# sequential grain colors
import colorsys
colorScale = (Vector3(colorsys.hsv_to_rgb(value*1.0/numThreads, 1, 1)) for value in range(0, numThreads))

#add spheres
for sd in range(0,numThreads-1):
	col = next(colorScale)
	if rank!=(sd+1): continue
	ids=[]
	_id = 0
	for i in range(N):#(numThreads-1) x N x M spheres, one thread is for master and will keep only the wall, others handle spheres
		for j in range(M):
			id = O.bodies.insertAtId(sphere((sd*N+i+j/30.,j,0),0.500,color=col),_id+(N*M*sd)) #a small shift in x-positions of the rows to break symmetry
			_id+=1
			ids.append(id)
	if rank is not None:# assigning subdomain!=0 in single thread would freeze the particles (Newton skips them)
		for id in ids: O.bodies[id].subdomain = sd+1

WALL_ID=O.bodies.insertAtId(box(center=(numThreads*N*0.5,-0.5,0),extents=(2*numThreads*N,0,2),fixed=True),1+(N*M*(numThreads-1)))

collider.verletDist = 0.25
newton.gravity=(0,-10,0) #else nothing would move
tsIdx=O.engines.index(timeStepper) #remove the automatic timestepper. Very important: we don't want subdomains to use many different timesteps...
O.engines=O.engines[0:tsIdx]+O.engines[tsIdx+1:]
O.dt=0.001 #this very small timestep will make it possible to run 2000 iter without merging
#O.dt=0.1*PWaveTimeStep() #very important, we don't want subdomains to use many different timesteps...


#########  RUN  ##########
def collectTiming():
	created = os.path.isfile("collect.dat")
	f=open('collect.dat','a')
	if not created: f.write("numThreads mpi omp Nspheres N M runtime \n")
	from yade import timing
	f.write(str(numThreads)+" "+str(os.getenv('OMPI_COMM_WORLD_SIZE'))+" "+os.getenv('OMP_NUM_THREADS')+" "+str(N*M*(numThreads-1))+" "+str(N)+" "+str(M)+" "+str(timing.runtime())+"\n")
	f.close()


if rank is None: #######  Single-core  ######
	O.timingEnabled=True
	O.run(NSTEPS,True)
	#print "num bodies:",len(O.bodies)
	from yade import timing
	timing.stats()
	collectTiming()
	print ("num. bodies:",len([b for b in O.bodies]),len(O.bodies))
	print ("Total force on floor=",O.forces.f(WALL_ID)[1])
else: #######  MPI  ######
	#import yade's mpi module
	from yade import mpy as mp
	# customize
	mp.MAX_RANK_OUTPUT=4
	mp.YADE_TIMING=True
	mp.DISTRIBUTED_INSERT=True
	mp.mpirun(1) #this is to eliminate initialization overhead in Cundall number and timings
	from yade import timing
	timing.reset()
	t1=time.time()
	mp.mpirun(NSTEPS)
	t2=time.time()
	mp.mprint("num. bodies:",len([b for b in O.bodies])," ",len(O.bodies))
	if rank==0:
		mp.mprint( "Total force on floor="+str(O.forces.f(WALL_ID)[1]))
		mp.mprint("CPU wall time for ",NSTEPS," iterations:",t2-t1,"; Cundall number = ",N*M*(numThreads-1)*NSTEPS/(t2-t1))
		collectTiming()
	else: mp.mprint( "Partial force on floor="+str(O.forces.f(WALL_ID)[1]))
	#mp.mergeScene()
	#if rank==0: O.save('mergedScene.yade')
	mp.MPI.Finalize()
exit()
