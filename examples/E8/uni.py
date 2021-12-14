"""
A simple example to generate uniaxial specimen using
Abaqus Python script feature

Youngung Jeong
youngung.jeong@gmail.com


abaqus cae -noGUI uni.py
"""

from abaqus import *
from caeModules import *
from abaqusConstants import *
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)
import sketch, part, regionToolset, os
import numpy as np

## import local site-packages to use my packages.

## __py_pack_dir__='/home/younguj/anaconda2/lib/python2.7/site-packages/'
__py_pack_dir__='c:/Users/user/AppData/Local/Programs/Python/Python37/lib/site-packages/'
import os
os.sys.path.append(__py_pack_dir__)
import abaquspy.sketches.E8, abaquspy.sketches.drawTool
import abaquspy.lib.sets
import abaquspy.lib.datums
tensileBar=abaquspy.sketches.E8.tensileBar
draw_arc=abaquspy.sketches.drawTool.draw_arc
draw_line=abaquspy.sketches.drawTool.draw_line
setNodeCoord=abaquspy.lib.sets.setNodeCoord

# Assign MidSpan using datum called datC: which is located
# in the center of specimen
def setTransSpan(myAssembly,edges,dat,name,tw):
    """
    Define a set attribute within myAssembly using
    the given datum <dat>. Select two neighboring edges that
    are aligned transverse to the axial direction of the specimen.
    Give the <name> to this set.

    Arguments
    ---------
    myAssembly
    edges
    dat
    name
    tw
    """
    coord1=np.array(dat.pointOn)
    coord2=np.array(dat.pointOn)
    coord1[1]=coord1[1]-tw/4.
    coord2[1]=coord2[1]+tw/4.
    C1=tuple(coord1);C2=tuple(coord2)
    myAssembly.Set(edges=(edges.findAt((C1,),(C2,))), name=name)

def returnDatumCoords(myPart,coords):
    ## below returns a feature
    f=myPart.DatumPointByCoordinate(coords=coords)
    ## convert feature to datum
    return myPart.datums[f.id]

## dimension in meter (default values are given in millimeters).
## thus a factor of 1e-3 is required.
#def main(Theta=0.,umatFN=None,myMatFunc=None,isub=False,iwait=False):
# """
# Arguments
# ---------
# Theta
# umatFN
# myMatFunc
# isub  = False
# iwait = False
# """
import logging
logging.basicConfig(filename='c:/users/user/repo/abaquspy/examples/E8/log.txt',
                    filemode='w',
                    datefmt='%H:%M:%S')
logging.info('** starts main')
logging.shutdown()

import sys


#def main(job_name='tmp',thick=2.,seedsize=5):

job_name='tmp'
thick=2.
seedsize=0.8  ##7884 # three layers
seedsize=1.0  ##3712 # two layers
#seedsize=2.0  ##0464 # single layer
#seedsize=3.0  ##0200 # single layer
#seedsize=5.0  ##0096 # single layer


Theta=0
label='%2.2i'%int(Theta)
label='%s_UMAT_None'%label

pl=57.
gw=12.5
tw=20.
tl=20.0
rd=12.5   ## Radius

xyCoords = tensileBar(pl=pl,gw=gw,tw=tw,tl=tl,rd=rd).T
xMax=np.max(xyCoords[:,0])
roundLength=(xMax-tl*2.-pl)/2.
rDC=returnDatumCoords # alias

myModel = mdb.Model(name='UT_%s'%label)
mySketch = myModel.ConstrainedSketch(name='E8_Sketch',sheetSize=thick)

## total length of specimen.
totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0])

for i in range(len(xyCoords)):
    mySketch.Spot(point=xyCoords[i])

draw_arc(mySketch, xyCoords,gw,0,0,1,'up')
draw_arc(mySketch, xyCoords,gw,6,5,6,'down')
draw_arc(mySketch, xyCoords,gw,8,8,9,'down')
draw_arc(mySketch, xyCoords,gw,14,13,14,'up')

draw_line(mySketch,xyCoords,1,5)
draw_line(mySketch,xyCoords,6,8)
draw_line(mySketch,xyCoords,9,13)
draw_line(mySketch,xyCoords,14,16)

###### ----------Generating specimen dimensions.
myPart = myModel.Part(name='E8', dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)

## Features.
if False:
    featShell=myPart.BaseShell(sketch=mySketch)
    myPart.features.changeKey(fromName=featShell.name,toName='myBaseShell')
else:
    featSolid=myPart.BaseSolidExtrude(sketch=mySketch,depth=2)
    myPart.features.changeKey(fromName=featSolid.name,toName='myBaseSolid')

L=totalLength
arcl=(L-pl-tl*2)/2. ## length of arc
gauge_bot_y = 0.5*tw-0.5*gw
gauge_top_y = 0.5*tw+0.5*gw
gby=gauge_bot_y
gty=gauge_top_y

faces=myPart.faces
cells=myPart.cells

rdc=returnDatumCoords # alias
## (a,b,c)
a=rdc(myPart,coords=(tl,tw,0))
b=rdc(myPart,coords=(tl,tw/2.,thick))
c=rdc(myPart,coords=(tl,0,0))
## (d,e,f)
d=rdc(myPart,coords=(tl+arcl,gty,0))
e=rdc(myPart,coords=(tl+arcl,tw/2.,thick))
f=rdc(myPart,coords=(tl+arcl,gby,0))
## (g,h,i)
g=rdc(myPart,coords=(tl+arcl+pl/4.,gty,0))
h=rdc(myPart,coords=(tl+arcl+pl/4.,tw/2.,thick))
i=rdc(myPart,coords=(tl+arcl+pl/4.,gby,0))
## (j,k,l)
j=rdc(myPart,coords=(tl+arcl+3*pl/4.,gty,0))
k=rdc(myPart,coords=(tl+arcl+3*pl/4.,tw/2.,thick))
l=rdc(myPart,coords=(tl+arcl+3*pl/4.,gby,0))
## (m,n,o)
m=rdc(myPart,coords=(tl+arcl+pl,gty,0))
n=rdc(myPart,coords=(tl+arcl+pl,tw/2.,thick))
o=rdc(myPart,coords=(tl+arcl+pl,gby,0))
## (p,q,r)
p=rdc(myPart,coords=(tl+arcl*2+pl,tw,0))
q=rdc(myPart,coords=(tl+arcl*2+pl,tw/2.,thick))
r=rdc(myPart,coords=(tl+arcl*2+pl,0,0))


# cell=myPart.cells
p1=myPart.DatumPlaneByThreePoints(point1=a,point2=b,point3=c)
p2=myPart.DatumPlaneByThreePoints(point1=d,point2=e,point3=f)
p3=myPart.DatumPlaneByThreePoints(point1=g,point2=h,point3=i)
p4=myPart.DatumPlaneByThreePoints(point1=j,point2=k,point3=l)
p5=myPart.DatumPlaneByThreePoints(point1=m,point2=n,point3=o)
p6=myPart.DatumPlaneByThreePoints(point1=p,point2=q,point3=r)

##-- set "entire specimen"
# c0=myPart.cells.findAt((0.,0,0))
#myPart.Set(cells=c0,name='Entire_Specimen')

## partitionning by planes
if True:
    c0=myPart.cells.findAt((0.0,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,a,b,c)
    c0=myPart.cells.findAt((L,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,d,e,f)
    c0=myPart.cells.findAt((L,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,g,h,i)
    c0=myPart.cells.findAt((L,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,j,k,l)
    c0=myPart.cells.findAt((L,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,m,n,o)
    c0=myPart.cells.findAt((L,0,0))
    myPart.PartitionCellByPlaneThreePoints(c0,p,q,r)

    myPart.seedPart(size=seedsize, deviationFactor=0.1, minSizeFactor=0.1)


    #elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)


    myPart.generateMesh()

    allcells=myPart.cells.findAt(
        ((0.01,tw/2.,0),),
        ((tl+0.01,tw/2.,0),),
        ((tl+arcl+0.01,tw/2.,0),),
        ((tl+arcl+pl/4.+0.01,tw/2.,0),),
        ((tl+arcl+3*pl/4.+0.01,tw/2.,0),),
        ((tl+arcl+pl+0.01,tw/2.,0),),
        ((L-0.01,tw/2.,0),),
    )


myMat = myModel.Material('myUMAT')
myMat.UserMaterial(mechanicalConstants=(0.,))
myMat.Depvar(n=7) ## Number of state variables


myMat_IFSTEEL = myModel.Material(name='ifsteel')
myMat_IFSTEEL.Plastic(table=((300.0, 0.0), (400.0, 0.1)))
myMat_IFSTEEL.Elastic(table=((200000.0, 0.3), ))
myMat_IFSTEEL.Density(table=((0.3,),))



a = myModel.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='myInstance', part=myPart, dependent=ON)
myInstance=a.instances['myInstance']


## generate faces to which boundary condition is applied
if True:
    face1=myInstance.faces.findAt(((0,tw/2.,thick/2.),))
    a.Set(faces=face1,name='leftEndFace')
    face1=myInstance.faces.findAt(((L,tw/2.,thick/2.),))
    a.Set(faces=face1,name='rightEndFace')
    verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    a.Set(vertices=verts1, name='Ori')

    # n1 = a.instances['myInstance'].nodes
    n1=myInstance.nodes.getClosest(coordinates=((tl+arcl+pl/4.,tw/2.,0),))
    n0=n1[0].label-1
    a.Set(nodes=myInstance.nodes[n0:n0+1],name='South')
    n1=myInstance.nodes.getClosest(coordinates=((tl+arcl+3*pl/4.,tw/2.,0),))
    n0=n1[0].label-1
    a.Set(nodes=myInstance.nodes[n0:n0+1],name='North')
    n1=myInstance.nodes.getClosest(coordinates=((tl+arcl+pl/2.,gty,0),))
    n0=n1[0].label-1
    a.Set(nodes=myInstance.nodes[n0:n0+1],name='West')
    n1=myInstance.nodes.getClosest(coordinates=((tl+arcl+pl/2.,gby,0),))
    n0=n1[0].label-1
    a.Set(nodes=myInstance.nodes[n0:n0+1],name='East')
    n1=myInstance.nodes.getClosest(coordinates=((tl+arcl+pl/2.,tw/2.,0),))
    n0=n1[0].label-1
    a.Set(nodes=myInstance.nodes[n0:n0+1],name='Center')







    # verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    # a.Set(vertices=verts1, name='North')
    # verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    # a.Set(vertices=verts1, name='South')
    # verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    # a.Set(vertices=verts1, name='East')
    # verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    # a.Set(vertices=verts1, name='West')
    # verts1 = myInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
    # a.Set(vertices=verts1, name='Center')




myModel.HomogeneousSolidSection(
    name='MySection',
    #material='myUMAT',
    material='ifsteel',
    thickness=None)

region=regionToolset.Region(cells=allcells)
myPart.SectionAssignment(region=region, sectionName='MySection', offset=0.0,
                    offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)

## Create a static general step
if True:
    myModel.StaticStep(
        name='Tension',previous='Initial',description='Uniaxial tension',
        timePeriod=1,adiabatic=OFF,maxNumInc=100,stabilization=None,
        timeIncrementationMethod=AUTOMATIC,initialInc=1e-1,minInc=1e-1,
        maxInc=1e-1,matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)
else:
    myModel.ImplicitDynamicsStep(
        name='Tension',previous='Initial',description='Uniaxial tension',
        timeIncrementationMethod=FIXED, initialInc=0.1,
        nohaf=OFF, noStop=OFF)


## history output
sets=['South','North','West','East','Center']
for s in sets:
    regionDef=a.sets[s]
    myModel.HistoryOutputRequest(
        name='H-Output-%s'%s,
        createStepName='Tension', variables=(
            'MISES', 'RF1', 'RF2', 'RF3', 'RM1',
            'RM2', 'RM3', 'RT', 'TF1', 'TF2', 'TF3', 'TM1', 'TM2', 'TM3'),
        region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)


myModel.DisplacementBC(
    name='BC-1',
    createStepName='Initial', region=a.sets['leftEndFace'], u1=0.0, u2=0, u3=0,
    ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF,
    distributionType=UNIFORM, fieldName='', localCsys=None)

region = a.sets['rightEndFace']
myModel.VelocityBC(
    name='BC-2', createStepName='Tension',
    region=region, v1=0.002, v2=0, v3=0, vr1=0.0, vr2=0.0, vr3=0.0,
    amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')

region = a.sets['Ori']
myModel.EncastreBC(
    name='BC-3', createStepName='Initial',
    region=region, localCsys=None)


##save cae.
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)
# mdb.saveAs(pathName='c:/Users/user/repo/abaquspy/examples/E8/ddum')

##
mdb.Job(name=job_name, model='UT_00_UMAT_None', description='', type=ANALYSIS,
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
    numGPUs=0)

mdb.jobs[job_name].writeInput(consistencyChecking=OFF)
mdb.jobs[job_name].submit(consistencyChecking=OFF)

def main_old(Theta=0.,umatFN=None,myMatFunc=None,isub=False,iwait=False):
    """
    Arguments
    ---------
    Theta
    umatFN
    myMatFunc
    isub  = False
    iwait = False
    """
    import logging
    logging.basicConfig(filename='c:/users/user/repo/abaquspy/examples/E8/log.txt',
                        filemode='w',
                        datefmt='%H:%M:%S')
    logging.info('** starts main')
    logging.shutdown()

    import sys
    label='%2.2i'%int(Theta)
    if type(umatFN)==type(None):
        label='%s_UMAT_None'%label
    else:
        label='%s_UMAT_%s'%(label,os.path.split(umatFN)[-1].split('.')[0])
    if type(myMatFunc)==type(None):
        label='%s_MatF_None'%label
    else:
        label='%s_MatF_%s'%(label,myMatFunc.__name__)

    mpa=1e6
    gpa=1e9
    mm = 1e-3
    pl=57.    * mm
    gw=12.5   * mm
    tw=20.    * mm
    tl=20.0   * mm
    rd=12.5   * mm ## Radius
    xyCoords = tensileBar(pl=pl,gw=gw,tw=tw,tl=tl,rd=rd).T
    sys.stdout.write('len(xyCoords):%i'%len(xyCoords))

    myModel = mdb.Model(name='UT_%s'%label)
    mySketch = myModel.ConstrainedSketch(name='E8_Sketch',sheetSize=1.0)

    ## total length of specimen.
    totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0])

    for i in range(len(xyCoords)):
        mySketch.Spot(point=xyCoords[i])

    draw_arc(mySketch, xyCoords,gw,0,0,1,'up')
    draw_arc(mySketch, xyCoords,gw,6,5,6,'down')
    draw_arc(mySketch, xyCoords,gw,8,8,9,'down')
    draw_arc(mySketch, xyCoords,gw,14,13,14,'up')
    #draw_arc(mySketch, xyCoords,gw,8,7,8,'down')
    #draw_arc(mySketch, xyCoords,gw,10,10,11,'down')
    #draw_arc(mySketch, xyCoords,gw,18,17,18,'up')

    draw_line(mySketch,xyCoords,1,5)
    draw_line(mySketch,xyCoords,6,8)
    draw_line(mySketch,xyCoords,9,13)
    draw_line(mySketch,xyCoords,14,16)

    ###### ----------Generating specimen dimensions.
    myPart = myModel.Part(name='E8', dimensionality=THREE_D,
                          type=DEFORMABLE_BODY)

    ## Features.
    if False:
        featShell=myPart.BaseShell(sketch=mySketch)
        myPart.features.changeKey(fromName=featShell.name,toName='myBaseShell')
    else:
        featSolid=myPart.BaseSolidExtrude(sketch=mySketch,depth=2)
        myPart.features.changeKey(fromName=featSolid.name,toName='myBaseSolid')

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

    rDC=returnDatumCoords # alias

    datOri=rDC(myPart,coords=(0,0,0))
    datO=rDC(myPart,coords=(totalLength      ,0,    0))
    datX=rDC(myPart,coords=(totalLength-tw/2.,0,    0))
    datY=rDC(myPart,coords=(totalLength      ,tw/2.,0))
    datC=rDC(myPart,coords=(totalLength/2.,tw/2.,0))
    datC_up=rDC(myPart,coords=(totalLength/2.,tw/2.+gw/2.,0))
    datC_down=rDC(myPart,coords=(totalLength/2.,tw/2.-gw/2.,0))
    datC_Lend=rDC(myPart,coords=(0,tw/2.,0))
    datC_Rend=rDC(myPart,coords=(totalLength,tw/2.,0))

    datC_Left=rDC(myPart,coords=(totalLength/2.-pl/2.,tw/2.,0))
    datC_LeftUp=rDC(myPart,coords=(totalLength/2.-pl/2.,tw/2.+gw/2.,0))
    datC_LeftDown=rDC(myPart,coords=(totalLength/2.-pl/2.,tw/2.-gw/2.,0))

    datC_Right=rDC(myPart,coords=(totalLength/2.+pl/2.,tw/2.,0))
    datC_RightUp=rDC(myPart,coords=(totalLength/2.+pl/2.,tw/2.+gw/2.,0))
    datC_RightDown=rDC(myPart,coords=(totalLength/2.+pl/2.,tw/2.-gw/2.,0))

    # Define Point Datums at the two edges of Grips.
    datPs=np.empty((2,3),dtype='object')
    xs=[tl,tl+pl+2.*roundLength]
    ys=[0,tw/2.,tw]
    for i in range(2):
        for j in range(3):
            datPs[i,j] = rDC(myPart=myPart,coords=(xs[i],ys[j],0))

    SysDefault = myPart.DatumCsysByDefault(CARTESIAN)
    cSysMat = abaquspy.lib.datums.returnCsymPlanarAngle(
        myPart,datC,radius=tw/3.,angle=Theta*np.pi/180.)

    ## Shell
    myPart.BaseShell(sketch=mySketch)
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

    ## Create My material
    if type(umatFN)==type(None):
        myMat = myModel.Material('IFsteel') ## modules:
        if type(myMatFunc)==type(None):
            abaquspy.mats.ifsteel.isoep(myMat)
        else:
            myMatFunc(myMat)
    else:
        myMat = myModel.Material('myUMAT')
        myMat.UserMaterial(mechanicalConstants=(200*gpa,0.3))
        myMat.Depvar(n=7) ## Number of state variables

    ## Create Shell Section!
    thickness=1e-3
    myShellSection=myModel.HomogeneousShellSection(
        name='SpecimenSection',preIntegrate=OFF, material=myMat.name,
        thicknessType=UNIFORM, thickness=thickness,thicknessField='',
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT,
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,
        integrationRule=SIMPSON, numIntPts=5)
    myShellSection.TransverseShearShell(
        k11=200.*gpa,k22=200.*gpa,k12=120.*gpa)

    ## Assign material orientation
    myPart.MaterialOrientation(localCsys=cSysMat,axis=AXIS_3)

    ## Assign section to plate
    region=regionToolset.Region(faces=myPart.faces[:])
    myPart.SectionAssignment(region=region,
                             sectionName=myShellSection.name)

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
    transPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datC_up,point2=datC_down)
    myPart.features.changeKey(fromName=transPlane.name,
                              toName='transPlane')
    #- axialPlane  -  along longitudinal direction of the specimen.
    axialPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datC_Lend,point2=datC_Rend)
    myPart.features.changeKey(fromName=axialPlane.name,
                              toName='axialPlane')
    #- pLeftPlane
    transPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datC_LeftDown,point2=datC_LeftUp)
    myPart.features.changeKey(fromName=transPlane.name,
                              toName='pLeftPlane')
    #- pRightPlane
    transPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datC_RightDown,point2=datC_RightUp)
    myPart.features.changeKey(fromName=transPlane.name,
                              toName='pRightPlane')
    #- gLeftPlane
    transPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datPs[0,0],point2=datPs[0,2])
    myPart.features.changeKey(fromName=transPlane.name,
                              toName='gLeftPlane')
    #- gRightPlane
    transPlane=myPart.PartitionFaceByShortestPath(
        faces=myPart.faces[:],point1=datPs[1,0],point2=datPs[1,2])
    myPart.features.changeKey(fromName=transPlane.name,
                              toName='gRightPlane')

    # MidSpan
    SpecimenNameInAssembly=list(myAssembly.instances.items())[0][0]
    edges=myAssembly.instances[SpecimenNameInAssembly].edges

    ## trans line sets.
    setTransSpan(myAssembly=myAssembly,edges=edges,dat=datC,
                 name='MidSpan',tw=tw)
    setTransSpan(myAssembly=myAssembly,edges=edges,dat=datC_Left,
                 name='ParallelLeft',tw=tw)
    setTransSpan(myAssembly=myAssembly,edges=edges,dat=datC_Right,
                 name='ParallelRight',tw=tw)
    setTransSpan(myAssembly=myAssembly,edges=edges,dat=datPs[0,1],
                 name='GripLeft',tw=tw)
    setTransSpan(myAssembly=myAssembly,edges=edges,dat=datPs[1,1],
                 name='GripRight',tw=tw)

    ## Create a static general step
    myModel.StaticStep(
        name='Tension',previous='Initial',description='Uniaxial tension',
        timePeriod=1,adiabatic=OFF,maxNumInc=100,stabilization=None,
        timeIncrementationMethod=AUTOMATIC,initialInc=1,minInc=1e-5,
        maxInc=1,matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)

    ## Define boundary conditions...
    epsRate=1.e-3 #0.001/sec
    delEpsMax=1.e-3 ## I want incremental step less than ...
    maxTimeInc=delEpsMax/epsRate
    minTimeInc=maxTimeInc/1000.
    # approximate gauge length:
    L0=0.95*pl
    vel=epsRate*L0 ## velocity

    ## total (engi) strain wanted: 0.2
    totalStrain = 0.10
    Lf=(1.+totalStrain)*L0
    totalDisplace=Lf-L0
    deltaTime=totalDisplace/vel ## total delta Time

    print('minInc:', minTimeInc)
    print('maxInc:', maxTimeInc)

    myModel.StaticStep(
        name='TensionContinue',previous='Tension',
        description='Uniaxial Tension',timePeriod=deltaTime,
        adiabatic=OFF,maxNumInc=1000,stabilization=None,
        timeIncrementationMethod=AUTOMATIC,initialInc=minTimeInc,
        minInc=minTimeInc,maxInc=maxTimeInc,matrixSolver=SOLVER_DEFAULT,
        extrapolation=DEFAULT)

    ## view
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Tension')

    ## Modify output request
    # Field output
    if type(umatFN)==type(None):
        myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('E','U','S','PE'))
    else:
        myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('E','U','S','UVARM'))

    # History output
    myModel.historyOutputRequests['H-Output-1'].setValues(
        variables=('E11',),region=myAssembly.sets['MidSpan'])

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON,
        predefinedFields=ON)

    ## Apply BC
    c0=datOri.pointOn
    vs=myInstance.vertices.findAt((c0,))
    myModel.EncastreBC(name='EncastreOri',createStepName='Tension',
                       region=regionToolset.Region(vertices=vs))
    myModel.XsymmBC(name='FixLeftEndX', createStepName='Tension',
                    region=myInstance.sets['leftEnd'])
    myModel.ZsymmBC(name='FixLeftEndZ', createStepName='Tension',
                    region=myInstance.sets['leftEnd'])
    myModel.ZsymmBC(name='FixRightEndZ',createStepName='Tension',
                    region=myInstance.sets['rightEnd'])

    # Velocity
    myModel.VelocityBC(name='StretchX', createStepName='Tension',
                       region=myInstance.sets['rightEnd'])
    myModel.boundaryConditions['StretchX'].setValuesInStep(
        stepName='TensionContinue',v1=vel*10,vr3=0.)

    ## Generate Mesh
    elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
    #elemType1 = mesh.ElemType(elemCode=S8R5,  elemLibrary=STANDARD))
    #elemType2 = mesh.ElemType(elemCode=STRI65,elemLibrary=STANDARD))

    f = myPart.faces
    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    pickedRegions =(faces, )
    myPart.setElementType(regions=pickedRegions,
                          elemTypes=(elemType1, elemType2))
    #myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
    myPart.seedPart(size=0.001, minSizeFactor=0.1) ## finer meshing
    myPart.generateMesh()

    ## Create Job
    jobName='UniE8_%s'%label
    mdb.Job(
        name=jobName,model=myModel.name,
        description='PythonScriptedUniaxialTensile_%s'%jobName)
    myJob=mdb.jobs[jobName]
    myAssembly.regenerate()
    setNodeCoord(myPart=myPart,dat=datC,name='Center',offset=1e-4)
    myAssembly.regenerate()
    mdb.saveAs(myModel.name)

    ## Flag to use a User Material subroutine
    if type(umatFN)!=type(None):
        print('User material has been specified.')
        myJob.setValues(userSubroutine=umatFN)
    if isub:
        ## submit the job
        myJob.submit(consistencyChecking=OFF)
        if iwait:
            ## myJob wait until completion?
            myJob.waitForCompletion()

    return myModel, myJob

#umatFN='/home/younguj/repo/abaqusPy/umats/epl/epl.f'
#umatFN='c:/Users/user/repo/abaquspy/umats/epl/epl.f'
#main_old(umatFN=umatFN,Theta=0.)
#main(job_name='tmp',seedsize=5.0)
