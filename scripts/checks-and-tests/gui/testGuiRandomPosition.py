#!/usr/bin/python
# -*- coding: utf-8 -*-

############################################################################################################################################
# FIXME NOTE: This is a verbatim copy of examples/ResetRandomPosition.py and O.run() removed from it.  Then the extra lines at the end
#             turn it into a GUI test. The only "human" work is to find proper values for O.dt and guiIterPeriod. Hence all GUI tests can be
#             largely automated with the example script taken and modified automatically. This auotmation should be done in the future.
############################################################################################################################################

from yade import pack,export,qt
import gts,os

def Plane(v1,v2,v3,v4):
	pts = [ [Vector3(v1),Vector3(v2),Vector3(v3),Vector3(v4)] ]
	return pack.sweptPolylines2gtsSurface(pts,capStart=True,capEnd=True)

# Parameters
tc=0.001# collision time 
en=0.3  # normal restitution coefficient
es=0.3  # tangential restitution coefficient
frictionAngle=radians(35)# 
density=2700
# facets material
facetMat=O.materials.append(ViscElMat(frictionAngle=frictionAngle,tc=tc,en=en,et=es)) 
# default spheres material
dfltSpheresMat=O.materials.append(ViscElMat(density=density,frictionAngle=frictionAngle,tc=tc,en=en,et=es))

O.dt=.2*tc # time step

Rs=0.02 # mean particle radius
Rf=0.01 # dispersion (Rs±Rf*Rs)
nSpheres=1000# number of particles

# Create geometry
pln=Plane( (-.5, -.5, 0), (.5, -.5, -.05), (.5, .5, 0), (-.5, .5, -.05) ); 
plnIds=O.bodies.append(pack.gtsSurface2Facets(pln,material=facetMat,color=(0,1,0)))

fct=Plane( (-.25, -.25, .5), (.25, -.25, .5), (.25, .25, .5), (-.25, .25, .5) ); 
fctIds=O.bodies.append(pack.gtsSurface2Facets(fct,material=facetMat,color=(1,0,0),noBound=True))

# Create spheres
sp=pack.SpherePack(); 
sp.makeCloud(Vector3(-.5, -.5, 0),Vector3(.5, .5, .2), Rs, Rf, int(nSpheres), False)
spheres=O.bodies.append([sphere(s[0],s[1],color=(0.929,0.412,0.412),material=dfltSpheresMat) for s in sp])

# Create engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_ViscElMat_ViscElMat_ViscElPhys()],
		[Law2_ScGeom_ViscElPhys_Basic()],
	),
	NewtonIntegrator(damping=0,gravity=[0,0,-9.81]),
	ResetRandomPosition(virtPeriod=0.01,factoryFacets=fctIds,velocity=(0,0,-2),subscribedBodies=spheres,point=(0,0,-.5),normal=(0,0,1),maxAttempts=100),
]

############################################################################################################################################
############################################################# test GUI #####################################################################
############################################################################################################################################
# here start changes of script simple-scene-energy-tracking.py, maybe later this duplication of code above can be removed.
# The code below, takes screenshot before and after each GUI action. And yade is hopefully not crashing in between.
# The test runs also on debug build, so anyway we should get a useful backtrace from gitlab-CI

from testGuiHelper import TestGUIHelper
# FIXME: it should deduce the name automatically, it's the end of the filename. See also testGui.sh
#        if you add a new file, you have to manually add it into scripts/checks-and-tests/gui/testGui.sh
#        or even better finally fix this FIXME, so that it works automatically
#        and just like scripts/checks-and-tests/checks/checkList.py is finding all the files to run.
# FIXME: it automatically should add itself to engines, now I have to do this by hand.
scr = TestGUIHelper("RandomPosition")
guiIterPeriod=5000
O.engines = O.engines+[PyRunner(iterPeriod=guiIterPeriod,command='scr.screenshotEngine()')]
O.dt = O.dt*0.1
O.run(guiIterPeriod*scr.getTestNum()+1)

