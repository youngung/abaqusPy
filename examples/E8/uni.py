"""
A simple example to generate uniaxial specimen using
Abaqus Python script feature

Youngung Jeong
youngung.jeong@gmail.com
"""

from abaqus import *
from caeModules import *
from abaqusConstants import *
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)
import sketch, part, regionToolset, os
import numpy as np

## import local site-packages to use my packages.
import os
os.sys.path.append('/home/younguj/anaconda2/lib/python2.7/site-packages/')
import abaquspy.sketches.E8, abaquspy.sketches.drawTool, abaquspy.lib.sets
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

# [SI unit system]
# length in meters]
# Pressure in Pa: N/m^2

## dimension in meter (default values are given in millimeters).
## thus a factor of 1e-3 is required.
def main(Theta=0.,umatFN=None,myMatFunc=None,isub=False,iwait=False):
    """
    Arguments
    ---------
    Theta
    umatFN
    myMatFunc
    isub  = False
    iwait = False
    """
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

    myModel = mdb.Model(name='UT_%s'%label)
    mySketch = myModel.ConstrainedSketch(name='E8_Sketch',sheetSize=1.0)

    ## total length of specimen.
    totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0])

    for i in xrange(len(xyCoords)):
        mySketch.Spot(point=xyCoords[i])

    draw_arc(mySketch, xyCoords,gw, 0,0,1,'up')
    draw_arc(mySketch, xyCoords,gw,8,7,8,'down')
    draw_arc(mySketch, xyCoords,gw,10,10,11,'down')
    draw_arc(mySketch, xyCoords,gw,18,17,18,'up')

    draw_line(mySketch,xyCoords,1,7)
    draw_line(mySketch,xyCoords,8,10)
    draw_line(mySketch,xyCoords,11,17)
    draw_line(mySketch,xyCoords,18,20)

    ###### ----------Generating specimen dimensions.
    myPart = myModel.Part(name='E8', dimensionality=THREE_D,
                          type=DEFORMABLE_BODY)

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
    for i in xrange(2):
        for j in xrange(3):
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
    SpecimenNameInAssembly=myAssembly.instances.items()[0][0]
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

    print 'minInc:', minTimeInc
    print 'maxInc:', maxTimeInc

    myModel.StaticStep(
        name='TensionContinue',previous='Tension',
        description='Uniaxial Tension',timePeriod=deltaTime,
        adiabatic=OFF,maxNumInc=500,stabilization=None,
        timeIncrementationMethod=AUTOMATIC,initialInc=minTimeInc,
        minInc=minTimeInc,maxInc=maxTimeInc,matrixSolver=SOLVER_DEFAULT,
        extrapolation=DEFAULT)

    ## view
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Tension')

    ## Modify output request
    # Field output
    myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('E','U','S','PE'))
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
    myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
    #myPart.seedPart(size=0.001, minSizeFactor=0.1) ## finer meshing
    myPart.generateMesh()

    ## Create Job
    jobName='UniE8_%s'%label
    mdb.Job(name=jobName,model=myModel.name,
            description='PythonScriptedUniaxialTensile_%s'%jobName)
    myJob=mdb.jobs[jobName]
    myAssembly.regenerate()
    setNodeCoord(myPart=myPart,dat=datC,name='Center',offset=1e-4)
    myAssembly.regenerate()
    mdb.saveAs(myModel.name)


    ## Flag to use a User Material subroutine
    if type(umatFN)!=type(None):
        print 'User material has been specified.'
        myJob.setValues(userSubroutine=umatFN)
    if isub:
        ## submit the job
        myJob.submit(consistencyChecking=OFF)
        if iwait:
            ## myJob wait until completion?
            myJob.waitForCompletion()

    return myModel, myJob


umatFN='/home/younguj/repo/abaqusPy/umats/epl/mises.f'
main(umatFN=umatFN,Theta=0.)
