"""
A simple example to generate one element model using
Abaqus Python script feature

Youngung Jeong
youngung.jeong@gmail.com
"""

## add the site-package of my own...
import os
import numpy as np

## Load my local packages by adding path_site_packge to system path
## I found this allows to load personal packages within Abaqus PDE.
## Yet, I don't know how reliable this method is.
path_site_packages='/home/younguj/anaconda2/lib/python2.7/site-packages/'
orig_path = os.sys.path[::]
os.sys.path.append(path_site_packages)
import abaquspy.mats.ifsteel
import abaquspy.lib.datums
rPDC=abaquspy.lib.datums.returnPointDatumCoords

from abaqus import *
from abaqusConstants import *
from caeModules import *
import sketch, part, regionToolset
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

###### ----------Generating specimen dimensions.
# [SI unit system]
# length in meters]
# Pressure in Pa: N/m^2
# dimension in meter (default values are given in millimeters).
# thus a factor of 1e-3 is required.
def square(length=1.e-3):
    """
    Return xy coordinates for a square
    with the given <length>

    Argument
    --------
    length
    """
    xyCoords=[[0,0],[1,0],[1,1],[0,1],[0,0]]
    xyCoords=np.array(xyCoords) * length
    return xyCoords

length=1.e-3
xyCoords=square(length=length)

def TensileOneElement(
    Theta=0.,isub=False,iwait=False,
    umatFN=None,myMatFunc=None,totalStrain=0.05):
    """
    Prep one shell element tests

    Arguments
    ---------
    Theta        = 0 : c.c.w in-plane rotation of material
                      from global axis 1
    isub        = False
    iwait       = False
    umatFN      = None
    myMatFunc   = None
    totalStrain = 0.05
    """
    ### Model declaration
    gpa=1.e9
    label='%2.2i'%int(Theta)

    if type(umatFN)==type(None):
        label='%s_UMAT_None'%label
    else:
        label='%s_UMAT_%s'%(label,os.path.split(umatFN)[-1].split('.')[0])

    if type(myMatFunc)==type(None):
        label='%s_MatF_None'%label
    else:
        label='%s_MatF_%s'%(label,myMatFunc.__name__)

    myModel = mdb.Model(name='OneEL_%s'%label)

    ### MySketch
    mySketch = myModel.ConstrainedSketch(name='Square_Sketch',
                                         sheetSize=length)
    ### Sketch details.
    ## Connecting the four corners.
    totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0])
    for i in xrange(len(xyCoords)): mySketch.Spot(point=xyCoords[i])
    for i in xrange(len(xyCoords)-1):
        p1=tuple(xyCoords[i])
        p2=tuple(xyCoords[i+1])
        mySketch.Line(point1=p1,point2=p2)

    ### Part module
    myPart = myModel.Part(name='oneElement',
                          dimensionality=THREE_D,
                          type=DEFORMABLE_BODY)
    ### Features
    ## Using the Sketch generate BasedShell.
    featShell=myPart.BaseShell(sketch=mySketch)
    myPart.features.changeKey(fromName=featShell.name,
                              toName='myBaseShell')
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

    ### Edge sets pertaining to Part
    edge0=myPart.edges.findAt((0,0.0001,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='leftEnd')
    edge0=myPart.edges.findAt((length,0.0001,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='rightEnd')

    ### Define various point Datums
    ## Coordinate datum for material orientation
    ## assuming that x//RD, y//TD and z//ND
    datOri= rPDC(myPart, coords=(0,                0,    0))
    dat10 = rPDC(myPart, coords=(length,           0,    0))
    datRC = rPDC(myPart, coords=(length,   length/2.,    0))
    dat11 = rPDC(myPart, coords=(length,      length,    0))
    datC  = rPDC(myPart, coords=(length/2.,length/2.,    0))

    ### Define Point Datums at the two edges of Grips.
    SysDefault = myPart.DatumCsysByDefault(CARTESIAN)
    cSysMat = abaquspy.lib.datums.returnCsymPlanarAngle(
        myPart,datC,radius=length/3.,angle=Theta*np.pi/180.)

    ### Shell
    myPart.BaseShell(sketch=mySketch)
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

    ### Create My material
    ## Use abaquspy.mats.ifsteel module to retrieve material
    ## characteristics of the IF steel


    if type(umatFN)==type(None):
        myMat = myModel.Material('IFsteel') ## modules:
        if type(myMatFunc)==type(None):
            abaquspy.mats.ifsteel.isoep(myMat)
        else:
            myMatFunc(myMat)
    else:
        myMat = myModel.Material('myUMAT')
        myMat.UserMaterial(mechanicalConstants=(200e9,0.3))
        myMat.Depvar(n=20) ## Number of state variables

    ### Create Shell Section!
    thickness = 1e-3 ## 1 mm thickness
    myShellSection=myModel.HomogeneousShellSection(
        name='SpecimenSection',preIntegrate=OFF, material=myMat.name,
        thicknessType=UNIFORM, thickness=thickness,thicknessField='',
        idealization=NO_IDEALIZATION,poissonDefinition=DEFAULT,
        thicknessModulus=None,temperature=GRADIENT,useDensity=OFF,
        integrationRule=SIMPSON, numIntPts=5)
    myShellSection.TransverseShearShell(
        k11=200.*gpa,k22=200.*gpa,k12=120.*gpa)


    ### Assign material orientation
    myPart.MaterialOrientation(localCsys=cSysMat,axis=AXIS_3)

    ### Assign section to plate
    region=regionToolset.Region(faces=myPart.faces[:])
    myPart.SectionAssignment(region=region,sectionName=myShellSection.name)

    ### retrieve assembly
    myAssembly = myModel.rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=myAssembly)
    ## Set cSys to myAssembly
    myAssembly.DatumCsysByDefault(CARTESIAN)
    ## Instance the part E8
    myAssembly.Instance(name='MySpecimen', part=myPart, dependent=ON)
    myInstance=myAssembly.instances['MySpecimen']

    ### Create a static general step
    #myModel.StaticStep(
    #    name='Tension',previous='Initial',description='Uniaxial tension',
    #    timePeriod=1,adiabatic=OFF,maxNumInc=100,stabilization=None,
    #    timeIncrementationMethod=AUTOMATIC,initialInc=1,minInc=1e-5,
    #    maxInc=1,matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)

    ### Define boundary conditions...
    epsRate=1e-3 #0.001/sec
    delEpsMax=1e-5 ## I want incremental step less than ...
    minTimeInc=delEpsMax/epsRate
    ## approximate gauge length:
    L0=1.*length   ## one element
    vel=epsRate*L0 ## velocity

    ### total (engi) strain wanted: 0.02
    # totalStrain = 0.20
    Lf=(1.+totalStrain)*L0
    totalDisplace=Lf-L0
    deltaTime=totalDisplace/vel ## total delta Time

    ##
    myModel.StaticStep(
        name='TensionContinue',previous='Initial',
        description='Uniaxial Tension',timePeriod=deltaTime,
        adiabatic=OFF,maxNumInc=2000,stabilization=None,
        timeIncrementationMethod=AUTOMATIC,initialInc=minTimeInc,
        minInc=minTimeInc,maxInc=minTimeInc*10.,
        matrixSolver=SOLVER_DEFAULT,extrapolation=DEFAULT)
    ## view
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        step='TensionContinue')

    ## Modify output request
    # Field output
    myModel.fieldOutputRequests['F-Output-1'].setValues(
        variables=('E','U','S','EE','PE'))
    # History output
    # myModel.historyOutputRequests['H-Output-1'].setValues(
    # variables=('E11',),region=myAssembly.sets['MidSpan'])

    ## Apply BC
    c0=datOri.pointOn
    vs=myInstance.vertices.findAt((c0,))
    ## Encastre
    myModel.EncastreBC(name='EncastreOri',createStepName='Initial',
                       region=regionToolset.Region(vertices=vs))
    ## Symmetric constraints
    myModel.XsymmBC(name='FixLeftEndX', createStepName='Initial',
                    region=myInstance.sets['leftEnd'])
    myModel.ZsymmBC(name='FixLeftEndZ', createStepName='Initial',
                    region=myInstance.sets['leftEnd'])
    myModel.ZsymmBC(name='FixRightEndZ',createStepName='Initial',
                    region=myInstance.sets['rightEnd'])
    ## Velocity
    myModel.VelocityBC(name='StretchX', createStepName='Initial',
                       region=myInstance.sets['rightEnd'])
    myModel.boundaryConditions['StretchX'].setValuesInStep(
        stepName='TensionContinue',v1=vel,vr3=0.)

    ## Generate Mesh
    elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
    faces=myPart.faces[:]
    pickedRegions =(faces, )
    myPart.setElementType(regions=pickedRegions, elemTypes=(
            elemType1, elemType2))
    myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
    myPart.generateMesh()

    ## additional procedure that requires 'meshed' information.
    myAssembly.regenerate()
    myAssembly.Set(elements=myInstance.elements[:],name='ORIGIN')
    myAssembly.regenerate()

    myModel.HistoryOutputRequest(
        name='StressStrain',rebar=EXCLUDE,
        createStepName='TensionContinue',variables=(
            'S11','E11','E22','PE11','PE22','EE11','EE22'),
        region=myAssembly.sets['ORIGIN'],sectionPoints=DEFAULT)
    myAssembly.regenerate()

    ## Create Job
    jobName='OneElement_%s'%label
    mdb.Job(name=jobName,model=myModel.name,
            description='PythonScriptedOneElement_%s'%jobName)
    mdb.saveAs(myModel.name)
    myJob = mdb.jobs[jobName]
    ## Multicore options
    # myJob.setValues(numCpus=4,numDomains=4)

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

## parametric usage of TensileOneElement
## below is application.
def runSingle(**kwargs):
    """
    Arguments
    ---------
    **kwargs    key-worded arguments passed to
                TensileOneElement
    """
    myModel, myJob = TensileOneElement(**kwargs)

def runTensions(nth, **kwargs):
    """
    Arguments
    ---------
    nth
    **kwargs
    """
    ths=np.linspace(0,90,nth)
    for i in xrange(nth):
        myModel, myJob = TensileOneElement(
            Theta=ths[i],**kwargs)

def runVarMats(**kwargs):
    """
    Arguments
    ---------
    **kwargs
    """
    myMatFuncs=[abaquspy.mats.ifsteel.isoe,
                abaquspy.mats.ifsteel.isoep]
    for imat in xrange(len(myMatFuncs)):
        runSingle(myMatFunc=myMatFuncs[imat],**kwargs)

if __name__=='main':
    print


## controlling job conditions
#umatFN=None
#umatFN='/home/younguj/repo/abaqusPy/umats/el/iso.f'
umatFN='/home/younguj/repo/abaqusPy/umats/epl/mises.f'

## Job testing methods
runSingle(umatFN=umatFN,iwait=True,isub=True,totalStrain=0.01)

## testing at various angles
#runTensions(nth=3,umatFN=umatFN,isub=False,iwait=False)

#runVarMats(umatFN=None,    isub=True)
#runVarMats(umatFN=umatFN,  isub=True)

os.sys.path=orig_path[::]
