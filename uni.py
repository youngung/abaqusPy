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
def drawE8_Sketch(
    pl=57,   ## Parallel length
    gw=12.5, ## Gauge width
    tw=20,   ## Grip width
    tl=50.0, ## Grip length
    rd=12.5, ## Radius
    ):
    """
    Draw E8 uniaxial tension specimen
    Units are in [mm]

    Arguments
    ---------
    pl : parallel length  (57   mm)
    gw : gauge width      (12.5 mm)
    tw : grip width       (20   mm)
    tl : grip length      (50.0 mm)
    rd : radius           (12.5 mm)
    """

    A  = tw/2. - gw/2.
    th = -90.-np.arccos((rd-A)/rd)*180./np.pi
    x  = rd * np.sin(th*np.pi/180.)

    ## Round...
    th0=-90
    th_delta=np.arccos((rd-A)/rd)*180./np.pi
    th1=th0+th_delta
    # ths=np.linspace(th0*np.pi/180.,th1*np.pi/180.)
    ths=[th0*np.pi/180.,th1*np.pi/180.]
    xs = rd*np.cos(ths)
    ys = rd*np.sin(ths)
    ## translate xs,ys
    xs = xs + (x-xs[-1])
    ys = ys + (-A+tw-ys[-1])
    xyRound=[xs.tolist(),ys.tolist()]
    ## parallel
    x0,y0=xs[0],ys[0]
    xyParallel = [[x0-0.5*pl,x0],[y0,y0]]
    ## Right grip
    XS=[x+tl,x+tl,x][::-1]
    YS=[-A+0.5*tw,-A+tw,-A+tw][::-1]
    xyRG=[XS,YS]

    x=xyParallel[0]+xyRound[0]+xyRG[0]
    y=xyParallel[1]+xyRound[1]+xyRG[1]

    xyCoords=np.array([x,y])

    # print xyCoords.shape

    ## translate the coordinate so that the center of gravity is (0,0)
    xyCoords[0]=xyCoords[0]-xyCoords[0][0]
    xyCoords[1]=xyCoords[1]-xyCoords[1][-1]
    # plot(xyCoords[0],xyCoords[1],'-')

    ## Apply 2-fold symmetry.
    sym0 =[[ 1,0],[0, 1]] ## Identical
    sym1 =[[-1,0],[0, 1]] ## Reflect y axis
    sym2 =[[ 1,0],[0,-1]] ## Reflect x axis
    sym3 =[[-1,0],[0,-1]] ## Reflect origin

    sym = np.array([sym0,sym2,sym3,sym1])
    # plot(xyCoords[0,0],xyCoords[1,0],'x')

    xyTot=[[],[]]
    for i in xrange(len(sym)):
        symOp = sym[i][:,:]# (2,2)
        temp = np.tensordot(symOp,xyCoords,axes=[1,0])
        if i==1 or i==3:
            temp[0]=temp[0][::-1]
            temp[1]=temp[1][::-1]
        elif i==0 or i==2:
            temp=temp[::]

        for j in xrange(len(temp[0])):
            xyTot[0].append(temp[0][j])
            xyTot[1].append(temp[1][j])

    xyTot=np.array(xyTot)


    x0=min(xyTot[0])
    y0=min(xyTot[1])+tw/2.

    xyTot[0] = xyTot[0] - x0
    xyTot[1] = xyTot[1] - y0

    return xyTot

# ---
def tensileBar(pl=57,   ## Parallel length
               gw=12.5, ## Gauge width
               tw=20,   ## Grip width
               tl=50.0, ## Grip length
               rd=12.5, ## Radius
               ):
    """
    Close end.
    """
    xy=drawE8_Sketch(pl,gw,tw,tl,rd)
    xy[1,:] = xy[1,:]+tw/2.
    xyt=xy.T
    x,y=xy
    XYS=[]
    for i in xrange(len(xyt)-1):
        x0,y0 = xyt[i]
        x1,y1 = xyt[i+1]
        if (x1-x0)==0 and (y1-y0)==0:
            pass
        else:
            XYS.append([x1,y1])

    XYS.append(XYS[0])
    xyt=np.array(XYS).T
    return xyt
## --

#[SI unit system]
# length in meters]
# Pressure in Pa: N/m^2

## dimension in meter (default values are given in millimeters).
## thus a factor of 1e-3 is required.
pl=57.    * 1e-3
gw=12.5   * 1e-3
tw=20.    * 1e-3
tl=50.0   * 1e-3
rd=12.5   * 1e-3 ## Radius
xyCoords = tensileBar(pl=pl,gw=gw,tw=tw,tl=tl,rd=rd).T

myModel = mdb.Model(name='UniaxialTension')
mySketch = myModel.ConstrainedSketch(name='E8_Sketch',sheetSize=1.0)

#s.ArcByCenterEnds(center=?, point1=?,point2?)
totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0]) ## total length of specimen.

for i in xrange(len(xyCoords)):
    mySketch.Spot(point=xyCoords[i])

g,v,d,c=mySketch.geometry, mySketch.vertices, mySketch.dimensions, mySketch.constraints

def draw_arc(iref,i0,i1,opt='up',):
    """
    arguments
    ---------
    iref: index of coordinates to which circel center is referenced.
    i0  : arc start index.
    i1  : arc end index
    opt : 'up' or 'down' depending on the reletive location of circle center to specimen.
    """
    xyArc=xyCoords[iref]

    if opt=='up':    centerArc = (xyArc[0],xyArc[1]+gw)
    elif opt=='down':centerArc = (xyArc[0],xyArc[1]-gw)
    else: raise IOError, 'Unexpected opt given'

    p1=(xyCoords[i0][0],xyCoords[i0][1])
    p2=(xyCoords[i1][0],xyCoords[i1][1])
    mySketch.ArcByCenterEnds(center=centerArc, point1=p1, point2=p2)

def draw_line(i0,i1):
    n=i1-i0
    for i in xrange(n):
        I0=i0+i
        I1=I0+1
        p1=(xyCoords[I0][0],xyCoords[I0][1])
        p2=(xyCoords[I1][0],xyCoords[I1][1])

        mySketch.Line(point1=p1,point2=p2)

draw_arc(0,0,1,'up')
draw_arc(8,7,8,'down')
draw_arc(10,10,11,'down')
draw_arc(18,17,18,'up')

draw_line(1,7)
draw_line(8,10)
draw_line(11,17)
draw_line(18,20)

###### ----------Generating specimen dimensions.
myPart = myModel.Part(name='E8', dimensionality=THREE_D,type=DEFORMABLE_BODY)

def findEdgeBetween(i0,i1):
    x0,y0=xyCoords[i0]
    x1,y1=xyCoords[i1]
    x=(x0+x1)/2.
    y=(y0+y1)/2.
    return myPart.edges.findAt(x,y,0.)


## Features.
featShell=myPart.BaseShell(sketch=mySketch)
myPart.features.changeKey(fromName=featShell.name,toName='myBaseShell')
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

edge0=myPart.edges.findAt((0,0.0001,0))
edge1=myPart.edges[edge0.index-1:edge0.index+1]
myPart.Set(edges=edge1,name='leftEnd')

xyt=np.array(xyCoords).T
xs=xyt[0]
xMax=np.max(xs)
edge0=myPart.edges.findAt((xMax,0.0001,0))
edge1=myPart.edges[edge0.index:edge0.index+2]
myPart.Set(edges=edge1,name='rightEnd')
# Find length of round
roundLength=(xMax-tl*2.-pl)/2.


## create datum for material orientation
# - create three datum points
# assuming that x//RD, y//TD and z//ND
def returnDatumCoords(coords):
    ## below returns a feature
    f=myPart.DatumPointByCoordinate(coords=coords)
    ## convert feature to datum
    return myPart.datums[f.id]
rDC=returnDatumCoords # alias

datOri=rDC(coords=(0,0,0))
datO=rDC(coords=(totalLength      ,0,    0))
datX=rDC(coords=(totalLength-tw/2.,0,    0))
datY=rDC(coords=(totalLength      ,tw/2.,0))
datC=rDC(coords=(totalLength/2., tw/2.,0))
datC_up=rDC(coords=(totalLength/2., tw/2.+gw/2.,0))
datC_down=rDC(coords=(totalLength/2., tw/2.-gw/2.,0))
datC_Lend=rDC(coords=(0,tw/2.,0))
datC_Rend=rDC(coords=(totalLength,tw/2.,0))

datC_Left=rDC(coords=(totalLength/2.-pl/2., tw/2.,0))
datC_LeftUp=rDC(coords=(totalLength/2.-pl/2., tw/2.+gw/2.,0))
datC_LeftDown=rDC(coords=(totalLength/2.-pl/2., tw/2.-gw/2.,0))

datC_Right=rDC(coords=(totalLength/2.+pl/2., tw/2.,0))
datC_RightUp=rDC(coords=(totalLength/2.+pl/2., tw/2.+gw/2.,0))
datC_RightDown=rDC(coords=(totalLength/2.+pl/2., tw/2.-gw/2.,0))

# Define Point Datums at the two edges of Grips.
datPs=np.empty((2,3),dtype='object')
xs=[tl,tl+pl+2.*roundLength]
ys=[0,tw/2.,tw]
for i in xrange(2):
    for j in xrange(3):
        datPs[i,j] = rDC(coords=(xs[i],ys[j],0))


#cSysDefault = myPart.DatumCsysByDefault(CARTESIAN)
cSysMat     = myPart.DatumCsysByThreePoints(CARTESIAN,origin=datC,point1=datC_Right,point2=datC_up,name='matCsys')

## Shell
myPart.BaseShell(sketch=mySketch)
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)


## Create My material
#Elasto-Plastic metal
#- Young's modulus for steel:
Young=200.*1e9#[Convert GPa to Pa - Pa is a SI unit: N/m^2]
# Yield strength: 400MPa

myMat=myModel.Material('Metal') ## moduls:
myMat.Elastic(table=((Young,0.30),))
myMat.Plastic(table=((400.E6, 0.0), (
    420.E6, 0.02), (500.E6, 0.2), (600.E6, 0.5)))

## Create Shell Section!
myShellSection=myModel.HomogeneousShellSection(
    name='SpecimenSection',
    preIntegrate=OFF, material='Metal', thicknessType=UNIFORM, thickness=1e-3,
    thicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT,
    thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,
    integrationRule=SIMPSON, numIntPts=5)

#myShellSection=myModel.HomogeneousSolidSection(
#    name='SpecimenSection',
#    material=myMat.name,thickness=1.e-3)


## Assign material orientation
myPart.MaterialOrientation(localCsys=cSysMat,axis=AXIS_3)

## Assign section to plate
region=regionToolset.Region(faces=myPart.faces[:])
myPart.SectionAssignment(region=region,sectionName=myShellSection.name)

## retrieve assembly
myAssembly = myModel.rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=myAssembly)
# Set cSys to myAssembly
myAssembly.DatumCsysByDefault(CARTESIAN)
# Instance the part E8
myAssembly.Instance(name='MySpecimen', part=myPart, dependent=ON)
myInstance=myAssembly.instances['MySpecimen']

## Partition the plate
#- transPlane  -  along transverse direction to the specimen axis
transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_up,point2=datC_down)
myPart.features.changeKey(fromName=transPlane.name,toName='transPlane')
#- axialPlane  -  along longitudinal direction of the specimen.
axialPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_Lend,point2=datC_Rend)
myPart.features.changeKey(fromName=axialPlane.name,toName='axialPlane')
#- pLeftPlane
transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_LeftDown,point2=datC_LeftUp)
myPart.features.changeKey(fromName=transPlane.name,toName='pLeftPlane')
#- pRightPlane
transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datC_RightDown,point2=datC_RightUp)
myPart.features.changeKey(fromName=transPlane.name,toName='pRightPlane')
#- gLeftPlane
transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datPs[0,0],point2=datPs[0,2])
myPart.features.changeKey(fromName=transPlane.name,toName='gLeftPlane')
#- gRightPlane
transPlane=myPart.PartitionFaceByShortestPath(faces=myPart.faces[:],point1=datPs[1,0],point2=datPs[1,2])
myPart.features.changeKey(fromName=transPlane.name,toName='gRightPlane')

# MidSpan
SpecimenNameInAssembly=myAssembly.instances.items()[0][0]
edges=myAssembly.instances[SpecimenNameInAssembly].edges

# Assign MidSpan using datum called datC: which is located in the center of specimen
def setTransSpan(dat,name):
    coord1=np.array(dat.pointOn)
    coord2=np.array(dat.pointOn)
    coord1[1]=coord1[1]-tw/4.
    coord2[1]=coord2[1]+tw/4.
    C1=tuple(coord1);C2=tuple(coord2)
    myAssembly.Set(edges=(edges.findAt((C1,),(C2,))), name=name)

setTransSpan(dat=datC,name='MidSpan')
setTransSpan(dat=datC_Left,name='ParallelLeft')
setTransSpan(dat=datC_Right,name='ParallelRight')

setTransSpan(dat=datPs[0,1],name='GripLeft')
setTransSpan(dat=datPs[1,1],name='GripRight')


## Create a static general step
myModel.StaticStep(name='Tension',previous='Initial',description='Uniaxial tension',
                   timePeriod=1,
                   adiabatic=OFF,maxNumInc=100,
                   stabilization=None,timeIncrementationMethod=AUTOMATIC,
                   initialInc=1,minInc=1e-5,maxInc=1,
                   matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)


## Define boundary conditions...
epsRate=1.e-3 #0.001/sec
delEpsMax=1.e-5 ## I want incremental step less than ...
minTimeInc=delEpsMax/epsRate
# approximate gauge length:
L0=0.95*pl
vel=epsRate*L0 ## velocity

## total (engi) strain wanted: 0.2
totalStrain = 0.2
Lf=(1.+totalStrain)*L0
totalDisplace=Lf-L0
deltaTime=totalDisplace/vel ## total delta Time

##

myModel.StaticStep(name='TensionContinue',previous='Tension',description='Uniaxial Tension',
                   timePeriod=deltaTime,
                   adiabatic=OFF,maxNumInc=1000,
                   stabilization=None,timeIncrementationMethod=AUTOMATIC,
                   initialInc=minTimeInc,minInc=minTimeInc,maxInc=minTimeInc*100.,
                   matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)
## view
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Tension')


## Modify output request
# Field output
myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('E','U','S','PE'))
# History output
myModel.historyOutputRequests['H-Output-1'].setValues(variables=('E11',),region=myAssembly.sets['MidSpan'])

session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON,
    predefinedFields=ON)


## Apply BC
c0=datOri.pointOn
vs=myInstance.vertices.findAt((c0,))
myModel.EncastreBC(name='EncastreOri',createStepName='Tension',region=regionToolset.Region(vertices=vs))

myModel.XsymmBC(name='FixLeftEndX', createStepName='Tension',region=myInstance.sets['leftEnd'])
myModel.ZsymmBC(name='FixLeftEndZ', createStepName='Tension',region=myInstance.sets['leftEnd'])
myModel.ZsymmBC(name='FixRightEndZ',createStepName='Tension',region=myInstance.sets['rightEnd'])

# Velocity
myModel.VelocityBC(name='StretchX', createStepName='Tension',region=myInstance.sets['rightEnd'])
myModel.boundaryConditions['StretchX'].setValuesInStep(
    stepName='TensionContinue',
    v1=vel*10,vr3=0.)


## Generate Mesh
elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
#elemType1 = mesh.ElemType(elemCode=S8R5)
#elemType2 = mesh.ElemType(elemCode=STRI65)

f = myPart.faces
faces = f.getSequenceFromMask(mask=('[#f ]', ), )
pickedRegions =(faces, )
myPart.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
#myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
myPart.seedPart(size=0.001, minSizeFactor=0.1) ## finer meshing
myPart.generateMesh()

## Create Job

mdb.Job(name='TensileE8',model=myModel.name,description='PythonScriptedUniaxialTensile')
myAssembly.regenerate()
mdb.saveAs(myModel.name)
