

# Possible executions of this script
# ./yadempi script.py #interactive will spawn 3 additional workers
# mpiexec -n 4 ./yadempi script.py #non interactive

NSTEPS=100 #turn it >0 to see time iterations, else only initilization TODO!HACK
#NSTEPS=50 #turn it >0 to see time iterations, else only initilization
N=50; M=50; #(columns, rows) per thread

import os
from yade import mpy as mp
numThreads = 4


#add spheres
young = 5e6
compFricDegree = 0.0
O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle = radians(compFricDegree), density= 2600, label='sphereMat'))
O.materials.append(FrictMat(young=young*100, poisson = 0.5, frictionAngle = compFricDegree, density =2600, label='wallMat'))

mn,mx=Vector3(0,0,0),Vector3(70,150,70)
pred = pack.inAlignedBox(mn,mx)
O.bodies.append(pack.regularHexa(pred,radius=2.80,gap=0, material='sphereMat'))


wallIds=aabbWalls([Vector3(-360,-1,-360),Vector3(360,360,360)],thickness=10.0, material='wallMat')
O.bodies.append(wallIds)

collider.verletDist = 2
newton.gravity=(0.05,-0.5,0.05) #else nothing would move
tsIdx=O.engines.index(timeStepper) #remove the automatic timestepper. Very important: we don't want subdomains to use many different timesteps...
O.engines=O.engines[0:tsIdx]+O.engines[tsIdx+1:]
O.dt=0.01

#########  RUN  ##########

def collectTiming():
	created = os.path.isfile("collect.dat")
	f=open('collect.dat','a')
	if not created: f.write("numThreads mpi omp Nspheres N M runtime \n")
	from yade import timing
	f.write(str(numThreads)+" "+str(os.getenv('OMPI_COMM_WORLD_SIZE'))+" "+os.getenv('OMP_NUM_THREADS')+" "+str(N*M*(numThreads-1))+" "+str(N)+" "+str(M)+" "+str(timing.runtime())+"\n")
	f.close()

# customize mpy

mp.VERBOSE_OUTPUT=False
mp.YADE_TIMING=False
mp.DOMAIN_DECOMPOSITION= True
mp.MERGE_W_INTERACTIONS=False
mp.ERASE_REMOTE_MASTER=False
mp.REALLOCATE_FREQUENCY=5
mp.mpirun(NSTEPS,4,True)

## single-thread vtk output from merged scene
#if mp.rank == 0:
	#from yade import export
	#v=export.VTKExporter("mpi3d")
	
#for k in range(300):
	#mp.mpirun(30,4,True)
	#if mp.rank == 0:
		#v.exportSpheres(what=dict(subdomain='b.subdomain'))
