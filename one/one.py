"""
A simple example to generate uniaxial specimen using
ABQUS Python scripting feature.


Youngung Jeong
youngung.jeong@gmail.com
"""

import os

from abaqus import *
from abaqusConstants import *
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

import sketch
import part
import regionToolset
import numpy as np

from caeModules import *


###### ----------Generating specimen dimensions.
#[SI unit system]
# length in meters]
# Pressure in Pa: N/m^2

## dimension in meter (default values are given in millimeters).
## thus a factor of 1e-3 is required.
def square(length=1.e-3):
    """
    """
    xyCoords=[[0,0],[1,0],[1,1],[0,1],[0,0]]
    xyCoords=np.array(xyCoords) * length
    return xyCoords

length=1.e-3
xyCoords=square(length=length)
    
myModel = mdb.Model(name='OneElementTest')
mySketch = myModel.ConstrainedSketch(name='Square_Sketch',sheetSize=length)

#s.ArcByCenterEnds(center=?, point1=?,point2?)
totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0]) ## total length of specimen.

for i in xrange(len(xyCoords)):
    mySketch.Spot(point=xyCoords[i])

g,v,d,c=mySketch.geometry, mySketch.vertices, mySketch.dimensions, mySketch.constraints

for i in xrange(len(xyCoords)-1):
    p1=tuple(xyCoords[i])
    p2=tuple(xyCoords[i+1])
    print p1,p2
    mySketch.Line(point1=p1,point2=p2)


###### ----------Generating specimen dimensions.
myPart = myModel.Part(name='oneElement', dimensionality=THREE_D,type=DEFORMABLE_BODY)

# def findEdgeBetween(i0,i1):
#     x0,y0=xyCoords[i0]
#     x1,y1=xyCoords[i1]
#     x=(x0+x1)/2.
#     y=(y0+y1)/2.
#     return myPart.edges.findAt(x,y,0.)


# ## Features.
featShell=myPart.BaseShell(sketch=mySketch)
myPart.features.changeKey(fromName=featShell.name,toName='myBaseShell')
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

edge0=myPart.edges.findAt((0,0.0001,0))
edges=myPart.edges[edge0.index:edge0.index+1]
myPart.Set(edges=edges,name='leftEnd')
edge0=myPart.edges.findAt((length,0.0001,0))
edges=myPart.edges[edge0.index:edge0.index+1]
myPart.Set(edges=edges,name='rightEnd')

# ## create datum for material orientation
# # - create three datum points
# # assuming that x//RD, y//TD and z//ND
def returnDatumCoords(coords):
    ## below returns a feature
    f=myPart.DatumPointByCoordinate(coords=coords)
    ## convert feature to datum
    return myPart.datums[f.id]
rDC=returnDatumCoords # alias

datOri= rDC(coords=(0,0,0))
dat10 = rDC(coords=(length,0,    0))
datRC = rDC(coords=(length,length/2.,    0))
dat11 = rDC(coords=(length,length,0))
datC=rDC(coords=(length/2., length/2.,0))
# datC_up=rDC(coords=(totalLength/2., tw/2.+gw/2.,0))
# datC_down=rDC(coords=(totalLength/2., tw/2.-gw/2.,0))
# datC_Lend=rDC(coords=(0,tw/2.,0))
# datC_Rend=rDC(coords=(totalLength,tw/2.,0))

# datC_Left=rDC(coords=(totalLength/2.-pl/2., tw/2.,0))
# datC_LeftUp=rDC(coords=(totalLength/2.-pl/2., tw/2.+gw/2.,0))
# datC_LeftDown=rDC(coords=(totalLength/2.-pl/2., tw/2.-gw/2.,0))

# datC_Right=rDC(coords=(totalLength/2.+pl/2., tw/2.,0))
# datC_RightUp=rDC(coords=(totalLength/2.+pl/2., tw/2.+gw/2.,0))
# datC_RightDown=rDC(coords=(totalLength/2.+pl/2., tw/2.-gw/2.,0))

# # Define Point Datums at the two edges of Grips.

SysDefault = myPart.DatumCsysByDefault(CARTESIAN)
cSysMat     = myPart.DatumCsysByThreePoints(CARTESIAN,origin=datC,
                                            point1=datRC,point2=dat11,name='matCsys')

## Shell
myPart.BaseShell(sketch=mySketch)
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)


# ## Create My material
# #Elasto-Plastic metal
# #- Young's modulus for steel:
Young=200.*1e9#[Convert GPa to Pa - Pa is a SI unit: N/m^2]
# # Yield strength: 400MPa 

myMat=myModel.Material('Metal') ## moduls:
myMat.Elastic(table=((Young,0.30),))
myMat.Plastic(table=((400.E6, 0.0), (420.E6, 0.02), (500.E6, 0.2), (600.E6, 0.5),))

# ## Create Shell Section!
thickness = 1e-3 ## 1 mm thickness
myShellSection=myModel.HomogeneousShellSection(
    name='SpecimenSection',preIntegrate=ON,
    material=myMat.name,thickness=1.e-3,
    poissonDefinition=DEFAULT, temperature=GRADIENT)

# ## Assign material orientation
myPart.MaterialOrientation(localCsys=cSysMat,axis=AXIS_3)

# ## Assign section to plate
region=regionToolset.Region(faces=myPart.faces[:])
myPart.SectionAssignment(region=region,sectionName=myShellSection.name)

# ## retrieve assembly
myAssembly = myModel.rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=myAssembly)
# Set cSys to myAssembly
myAssembly.DatumCsysByDefault(CARTESIAN)
# Instance the part E8
myAssembly.Instance(name='MySpecimen', part=myPart, dependent=ON)
myInstance=myAssembly.instances['MySpecimen']

# ## Partition the plate
# #- transPlane  -  along transverse direction to the specimen axis
# transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_up,point2=datC_down)
# myPart.features.changeKey(fromName=transPlane.name,toName='transPlane')
# #- axialPlane  -  along longitudinal direction of the specimen.
# axialPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_Lend,point2=datC_Rend)
# myPart.features.changeKey(fromName=axialPlane.name,toName='axialPlane')
# #- pLeftPlane
# transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_LeftDown,point2=datC_LeftUp)
# myPart.features.changeKey(fromName=transPlane.name,toName='pLeftPlane')
# #- pRightPlane
# transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_RightDown,point2=datC_RightUp)
# myPart.features.changeKey(fromName=transPlane.name,toName='pRightPlane')
# #- gLeftPlane
# transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datPs[0,0],point2=datPs[0,2])
# myPart.features.changeKey(fromName=transPlane.name,toName='gLeftPlane')
# #- gRightPlane
# transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datPs[1,0],point2=datPs[1,2])
# myPart.features.changeKey(fromName=transPlane.name,toName='gRightPlane')

# # MidSpan
# SpecimenNameInAssembly=myAssembly.instances.items()[0][0]
# edges=myAssembly.instances[SpecimenNameInAssembly].edges

# # Assign MidSpan using datum called datC: which is located in the center of specimen
# def setTransSpan(dat,name):
#     coord1=np.array(dat.pointOn)
#     coord2=np.array(dat.pointOn)
#     coord1[1]=coord1[1]-tw/4.
#     coord2[1]=coord2[1]+tw/4.
#     C1=tuple(coord1);C2=tuple(coord2)
#     myAssembly.Set(edges=(edges.findAt((C1,),(C2,))), name=name)

# setTransSpan(dat=datC,name='MidSpan')
# setTransSpan(dat=datC_Left,name='ParallelLeft')
# setTransSpan(dat=datC_Right,name='ParallelRight')

# setTransSpan(dat=datPs[0,1],name='GripLeft')
# setTransSpan(dat=datPs[1,1],name='GripRight')


# ## Create a static general step
myModel.StaticStep(name='Tension',previous='Initial',description='Uniaxial tension',
                   timePeriod=1,
                   adiabatic=OFF,maxNumInc=100,
                   stabilization=None,timeIncrementationMethod=AUTOMATIC,
                   initialInc=1,minInc=1e-5,maxInc=1,
                   matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)


# ## Define boundary conditions...
epsRate=1e-3 #0.001/sec
delEpsMax=1e-4 ## I want incremental step less than ...
minTimeInc=delEpsMax/epsRate
# approximate gauge length:
#L0=0.95*pl
L0=1.*length ## one element
vel=epsRate*L0 ## velocity

# ## total (engi) strain wanted: 0.2
totalStrain = 0.5
Lf=(1.+totalStrain)*L0
totalDisplace=Lf-L0
deltaTime=totalDisplace/vel ## total delta Time

##

myModel.StaticStep(name='TensionContinue',previous='Tension',description='Uniaxial Tension',
                   timePeriod=deltaTime,
                   adiabatic=OFF,maxNumInc=2000,
                   stabilization=None,timeIncrementationMethod=AUTOMATIC,
                   initialInc=minTimeInc,minInc=minTimeInc,maxInc=minTimeInc*10.,
                   matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)
## view
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Tension')


## Modify output request
# Field output
myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('E','U','S'))
# History output
#myModel.historyOutputRequests['H-Output-1'].setValues(variables=('E11',),region=myAssembly.sets['MidSpan'])

# session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON,
#     predefinedFields=ON)


# ## Apply BC
c0=datOri.pointOn
vs=myInstance.vertices.findAt((c0,))
myModel.EncastreBC(name='EncastreOri',createStepName='Tension',region=regionToolset.Region(vertices=vs))

myModel.XsymmBC(name='FixLeftEndX', createStepName='Tension',region=myInstance.sets['leftEnd'])
myModel.ZsymmBC(name='FixLeftEndZ', createStepName='Tension',region=myInstance.sets['leftEnd'])
myModel.ZsymmBC(name='FixRightEndZ',createStepName='Tension',region=myInstance.sets['rightEnd'])

# # Velocity
myModel.VelocityBC(name='StretchX', createStepName='Tension',region=myInstance.sets['rightEnd'])
myModel.boundaryConditions['StretchX'].setValuesInStep(
    stepName='TensionContinue',
    v1=vel,vr3=0.)


# ## Generate Mesh
elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
f = myPart.faces
faces = f.getSequenceFromMask(mask=('[#f ]', ), )
pickedRegions =(faces, )
myPart.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
#myPart.seedPart(size=0.001, minSizeFactor=0.1) ## finer meshing
myPart.generateMesh()

## Create Job

mdb.Job(name='OneElement',model=myModel.name,description='PythonScriptedOneElement')
myAssembly.regenerate()
mdb.saveAs(myModel.name)
