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

def set_EncastreOri(
    myModel=None,myInstance=None,dat=None,
    createStepName='Initial',name='None'):
    """
    Arguments
    ---------
    myModel        = Abaqus Model
    dat            = point datum
    createStepName = 'initial'
    """
    c0=dat.pointOn
    vs=myInstance.vertices.findAt((c0,))
    ## Encastre
    myModel.EncastreBC(name=name,createStepName=createStepName,
                       region=regionToolset.Region(vertices=vs))

def add_user_orien(myPart,region):
    """
    Add user orietation to the defined <region> pertaining to <myPart>
    """
    myPart.MaterialOrientation(
        region=region, orientationType=USER,
        additionalRotationType=ROTATION_NONE,
        localCsys=None, fieldName='')


length=1.e-3
xyCoords=square(length=length)

def TensileOneElement(
    Theta=0.,isub=False,iwait=False,
    umatFN=None,myMatFunc=None,totalStrain=0.05,iload=0):
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
    iload       = 0:  uniaxial along X
                  1:  pure shear
                  2:  simple shear
                  3:  uniaxial along Y
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

    ## Model name
    myModelName='OneEL_%s'%label

    myModel = mdb.Model(name=myModelName)

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
    myPartName='oneElement'
    myPart = myModel.Part(name=myPartName,
                          dimensionality=THREE_D,
                          type=DEFORMABLE_BODY)
    ### Features
    ## Using the Sketch generate BasedShell.
    featShell=myPart.BaseShell(sketch=mySketch)
    myPart.features.changeKey(fromName=featShell.name,
                              toName='myBaseShell')
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

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
    dat01  = rPDC(myPart, coords=(0,length,    0))

    ### Edge sets pertaining to Part
    edge0=myPart.edges.findAt((0,0.0001,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='leftEnd')
    edge0=myPart.edges.findAt((length,0.0001,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='rightEnd')
    edge0=myPart.edges.findAt((0.0001,length,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='upEnd')
    edge0=myPart.edges.findAt((0.0001,0,0))
    edges=myPart.edges[edge0.index:edge0.index+1]
    myPart.Set(edges=edges,name='downEnd')
    vts = myPart.vertices.findAt(dat11.pointOn)
    vts = myPart.vertices[vts.index:vts.index+1]
    myPart.Set(vertices=vts,name='v11')
    vts = myPart.vertices.findAt(dat01.pointOn)
    vts = myPart.vertices[vts.index:vts.index+1]
    myPart.Set(vertices=vts,name='v01')

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
        myMat.UserOutputVariables(n=20)


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

    ### Assign material orientation
    region = regionToolset.Region(faces=myPart.faces[:])
    orientation =  mdb.models[myModelName].parts[myPartName].datums[cSysMat.id]
    mdb.models[myModelName].parts[myPartName].MaterialOrientation(
        region=region, orientationType=SYSTEM, axis=AXIS_3,
        localCsys=orientation,fieldName='',
        additionalRotationType=ROTATION_NONE, angle=0.0,
        additionalRotationField='')

    ### Define boundary conditions...
    epsRate=1e-3 #0.001/sec
    delEpsMax=1e-5 ## I want incremental step less than ...
    minTimeInc=delEpsMax/epsRate
    ## approximate gauge length:
    L0=1.*length   ## one element
    vel=epsRate*L0 ## velocity

    ### total (engi) strain wanted: 0.02
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
    if type(umatFN)==type(None):
        myModel.fieldOutputRequests['F-Output-1'].setValues(
            variables=('E','U','S','EE','PE'),directions=OFF)
    else:
        myModel.fieldOutputRequests['F-Output-1'].setValues(
            variables=('E','U','S','UVARM'),directions=OFF)

    ## Apply BC
    ## Encastre
    if iload==0: ## uniaxial loading along global direction X
        set_EncastreOri(myModel,myInstance,datOri,'Initial','fix00')
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
    elif iload==1:
        set_EncastreOri(myModel,myInstance,datOri,'Initial','fix00')
        set_EncastreOri(myModel,myInstance,dat10, 'Initial','fix10')
        ## Velocity v11
        myModel.VelocityBC(name='StretchX1', createStepName='Initial',
                           region=myInstance.sets['v11'])
        myModel.boundaryConditions['StretchX1'].setValuesInStep(
            stepName='TensionContinue',v1=vel,vr3=0.)
        ## Velocity v01
        myModel.VelocityBC(name='StretchX2', createStepName='Initial',
                           region=myInstance.sets['v01'])
        myModel.boundaryConditions['StretchX2'].setValuesInStep(
            stepName='TensionContinue',v1=vel,vr3=0.)
    elif iload==2:
        set_EncastreOri(myModel,myInstance,datOri,'Initial','fix00')
        set_EncastreOri(myModel,myInstance,dat10, 'Initial','fix10')
        ## Velocity v11
        myModel.VelocityBC(name='StretchX1', createStepName='Initial',
                           region=myInstance.sets['v11'])
        myModel.boundaryConditions['StretchX1'].setValuesInStep(
            stepName='TensionContinue',v1=vel,v2=0.,vr3=0.)
        ## Velocity v01
        myModel.VelocityBC(name='StretchX2', createStepName='Initial',
                           region=myInstance.sets['v01'])
        myModel.boundaryConditions['StretchX2'].setValuesInStep(
            stepName='TensionContinue',v1=vel,v2=0.,vr3=0.)
    elif iload==3: ## uniaxial tension along global direction Y
        set_EncastreOri(myModel,myInstance,datOri,'Initial','fix00')
        myModel.YsymmBC(name='FixDownEndY', createStepName='Initial',
                        region=myInstance.sets['downEnd'])
        myModel.ZsymmBC(name='FixDownEndZ', createStepName='Initial',
                        region=myInstance.sets['downEnd'])
        myModel.ZsymmBC(name='FixUpEndZ', createStepName='Initial',
                        region=myInstance.sets['upEnd'])
        ## Velocity
        myModel.VelocityBC(name='StretchY', createStepName='Initial',
                           region=myInstance.sets['upEnd'])
        myModel.boundaryConditions['StretchY'].setValuesInStep(
            stepName='TensionContinue',v2=vel,vr3=0.)



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
            'S11','E11','UVARM'),region=myAssembly.sets['ORIGIN'],
        sectionPoints=DEFAULT)
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
    import argparse
    parser.add_argument('-iopt',type=str,help='Options')
    parser.add_argument(
        '-umat',type=str,help='Options',
        default='/home/younguj/repo/abaqusPy/umats/epl/epl.f')
    parser.add_argument(
        '-iumat',action='store_true',
        help='whether or not use the umat')
    parser.add_argument(
        '-iwait',action='store_true',help='iwait')
    parser.add_argument(
        '-isub',action='store_true',help='isub')
    parser.add_argument(
        '-mxe',type=float,help='Maximum strain to be imposed',
        default=0.01)
    parser.add_argument(
        '-nth',type=float,help='Number of thetas',
        default=3)
    parser.add_argument(
        '-iload',type=float,help='Type of loading condition (0: unix, 1: ps, 2: ss, 3:uniy',
        default=3)
    args = parser.parse_args()

    if args.iopt==0: ## run at single material direction
        if args.iumat: umatFN=args.umat
        else: umatFN=None
        runSingle(umatFN=umatFN,
                  iwait=args.iwait,
                  isub=args.isub,
                  totalStrain=args.mxe)
    elif args.iopt==1: ## run at multiple number of directions
        if args.iumat: umatFN=args.umat
        else: umatFN=None
        runTensions(umatFN=umatFN,
                    iwait=args.iwait,
                    isub=args.isub,
                    nth=args.nth,
                    totalStrain=args.mxe)


### controlling job conditions
##umatFN=None
##umatFN='/home/younguj/repo/abaqusPy/umats/el/iso.f'
umatFN='/home/younguj/repo/abaqusPy/umats/epl/epl.f'

### Job testing methods
# runSingle(umatFN=umatFN,iwait=True,isub=True,totalStrain=0.01)

### testing at various angles
runTensions(nth=3,umatFN=umatFN,isub=True,iwait=True,totalStrain=0.01)

##runVarMats(umatFN=None,    isub=True)
##runVarMats(umatFN=umatFN,  isub=True)
os.sys.path=orig_path[::]
