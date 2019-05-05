# -*- coding: utf-8 -*-
# https://gitlab.com/yade-dev/trunk/issues/7
# https://answers.launchpad.net/yade/+question/220785
# https://bugs.launchpad.net/yade/+bug/1112763

from yade.utils import *

DbigSphere	= 1.
radius_bigSphere= DbigSphere/2
radius_sph	= DbigSphere/10

# Periodic cell
cellSizes=[DbigSphere*1. , DbigSphere*3. , DbigSphere*3. ]
# NOTE, these settings do not trigger the bug and work fine:
#cellSizes[0]= DbigSphere*1.000000000000001 # ok - Anything bigger  than that works fine
#cellSizes[0]= DbigSphere*0.605             # ok - Anything smaller than that works fine

failCount = 0
runCount  = 0

def printOK(sign,ordering,axis,initialSortAxis):
	pass
	#axis = (2-axis)%3 # bacause the movement is along Z axis, while axis numbering referes to X. Better to print the axis which causes problems, than axis used in the loop numbering.
	#print("\033[92m OK, sign=",sign," ordering=",ordering," axis=",axis," initialSortAxis=",initialSortAxis,"\033[0m")

def noteFail(sign,ordering,axis,initialSortAxis):
	global failCount
	failCount += 1
	axis = (2-axis)%3 # bacause the movement is along Z axis, while axis numbering referes to X. Better to print the axis which causes problems, than axis used in the loop numbering.
	print("\033[91m failure, sign=",sign," ordering=",ordering," axis=",axis," initialSortAxis=",initialSortAxis,"\033[0m")


for initialSortAxis in [0,1,2]:
	for axis in [0,1,2]:
		for sign in [-1.0,1.0]:
			for ordering in [True,False]:
				O.periodic  = True
				O.cell.setBox(Vector3(cellSizes[(0+axis)%3],cellSizes[(1+axis)%3],cellSizes[(2+axis)%3]))
				# Big blue sphere
				centerPos        =        [cellSizes[0]/2.,cellSizes[1]/2.,cellSizes[2]*(0.7 if sign>0 else 0.3)]
				center_bigSphere = Vector3(centerPos[(0+axis)%3],centerPos[(1+axis)%3],centerPos[(2+axis)%3]) # Vector3(cellSizes[0]/2.,cellSizes[1]/2.,cellSizes[2]*(0.7 if sign>0 else 0.3))
				bigSphere        = utils.sphere(center=center_bigSphere, radius=radius_bigSphere, color=(0,0,1))
				# Red small sphere
				centerSmall      =        [centerPos[0],	centerPos[1]+1./5.*DbigSphere,	centerPos[2]+sign*(-DbigSphere/2.-radius_sph)]
				center_sph	 = Vector3(centerSmall[(0+axis)%3],centerSmall[(1+axis)%3],centerSmall[(2+axis)%3]) # Vector3(center_bigSphere[0],center_bigSphere[1]+1./5.*DbigSphere,center_bigSphere[2]+sign*(-DbigSphere/2.-radius_sph))
				sphere		 = utils.sphere(center=center_sph, radius=radius_sph, color=(1,0,0))
				forceAxis        = [0,0,-sign*100]

				if(ordering==True):
					O.bodies.append((sphere,bigSphere))        # ok
				else:
					O.bodies.append((bigSphere,sphere))        # does NOT work

				O.engines=[
					ForceResetter(),
					InsertionSortCollider ([Bo1_Sphere_Aabb()], allowBiggerThanPeriod=True),
					InteractionLoop ([Ig2_Sphere_Sphere_ScGeom()],[Ip2_FrictMat_FrictMat_MindlinPhys(),],[Law2_ScGeom_MindlinPhys_Mindlin()]),
					ForceEngine (force=(forceAxis[(0+axis)%3],forceAxis[(1+axis)%3],forceAxis[(2+axis)%3]),ids=[bigSphere.id]),
					NewtonIntegrator (damping=0.5),
				]
				O.dt=0.2*utils.PWaveTimeStep()
				testedCollider=typedEngine("InsertionSortCollider")
				testedCollider.sortAxis=initialSortAxis

#
# uncomment line below if you run this script manually: yade ./scripts/checks-and-tests/checks/checkColliderByGiulia.py
#
				#input("\033[33mTo see the bug, uncomment line 388 in pkg/common/InsertionSortCollider.cpp and recompile.\033[0m\nYou can press ▶ (the start button), but to see something use x10 smaller timestep.\nWhen you finish examination press enter to continue this script. BTW: it will fail if you changed the timestep.\n")
				# at O.iter==2 is appears first time and i.isReal == False, at O.iter==4240, i.isReal==True and Overlap 1.1326899061803175e-05
				O.run(2,True)
#				print("dumpBounds():",testedCollider.dumpBounds()[(2-axis)%3])
				try:
					i=O.interactions[bigSphere.id,sphere.id]
					if((i.isReal == False) and (O.iter==2)):
						O.run(4240-3,True)
						wasReal = i.isReal
						O.run(1,True)
						if((wasReal==False) and (i.isReal == True) and (O.iter==4240) and (abs(O.interactions[bigSphere.id,sphere.id].geom.penetrationDepth - 1.1326899061803175e-05) < 1e-15 )):
							printOK(sign,ordering,axis,initialSortAxis)
							#print('iter first interaction\tStep',O.iter,'\tOverlap',O.interactions[bigSphere.id,sphere.id].geom.penetrationDepth)
						else:
							noteFail(sign,ordering,axis,initialSortAxis)
					else:
						noteFail(sign,ordering,axis,initialSortAxis)
				except Exception as e:
					noteFail(sign,ordering,axis,initialSortAxis)

				O.reset()
				runCount += 1

if failCount != 0:
	# https://misc.flogisoft.com/bash/tip_colors_and_formatting
	print("\033[93m checkColliderByGiulia, runCount=",runCount," failCount=",failCount,".\033[0m")
	resultStatus+=1

