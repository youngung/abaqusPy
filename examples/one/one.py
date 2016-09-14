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

##
Theta = 45. ## [in degree] cSym is ccw from global (lab) axis X


###### ----------Generating specimen dimensions.
# [SI unit system]
# length in meters]
# Pressure in Pa: N/m^2
# dimension in meter (default values are given in millimeters).
# thus a factor of 1e-3 is required.
def square(length=1.e-3):
    """
    """
    xyCoords=[[0,0],[1,0],[1,1],[0,1],[0,0]]
    xyCoords=np.array(xyCoords) * length
    return xyCoords
length=1.e-3
xyCoords=square(length=length)

#### Model declaration
myModel = mdb.Model(name='OneElementTest')

#### MySketch
mySketch = myModel.ConstrainedSketch(name='Square_Sketch',sheetSize=length)
### Sketch details.
## Connecting the four corners.
totalLength=np.max(xyCoords[:,0])-np.min(xyCoords[:,0]) ## total length of specimen.
for i in xrange(len(xyCoords)): mySketch.Spot(point=xyCoords[i])
#g,v,d,c=mySketch.geometry, mySketch.vertices, mySketch.dimensions, mySketch.constraints
for i in xrange(len(xyCoords)-1):
    p1=tuple(xyCoords[i])
    p2=tuple(xyCoords[i+1])
    print p1,p2
    mySketch.Line(point1=p1,point2=p2)


#### Part module
myPart = myModel.Part(name='oneElement',
                      dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)
### Features
# Using the Sketch generate BasedShell.
featShell=myPart.BaseShell(sketch=mySketch)
myPart.features.changeKey(fromName=featShell.name,toName='myBaseShell')
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

## Edge sets pertaining to Part
edge0=myPart.edges.findAt((0,0.0001,0))
edges=myPart.edges[edge0.index:edge0.index+1]
myPart.Set(edges=edges,name='leftEnd')
edge0=myPart.edges.findAt((length,0.0001,0))
edges=myPart.edges[edge0.index:edge0.index+1]
myPart.Set(edges=edges,name='rightEnd')

## Define various point Datums
# Coordinate datum for material orientation
# assuming that x//RD, y//TD and z//ND
datOri= rPDC(myPart, coords=(0,                0,    0))
dat10 = rPDC(myPart, coords=(length,           0,    0))
datRC = rPDC(myPart, coords=(length,   length/2.,    0))
dat11 = rPDC(myPart, coords=(length,      length,    0))
datC  = rPDC(myPart, coords=(length/2.,length/2.,    0))


## Define Point Datums at the two edges of Grips.
SysDefault = myPart.DatumCsysByDefault(CARTESIAN)
cSysMat = abaquspy.lib.datums.returnCsymPlanarAngle(
    myPart,datC,radius=length/3.,angle=Theta*np.pi/180.)
# cSysMat    = myPart.DatumCsysByThreePoints(
#     CARTESIAN,origin=datC,
#     point1=datRC,point2=dat11,name='matCsys')

## Shell
myPart.BaseShell(sketch=mySketch)
session.viewports['Viewport: 1'].setValues(displayedObject=myPart)

#### Create My material
## Use abaquspy.mats.ifsteel module to retrieve material characteristics
## of the IF steel
myMat = myModel.Material('IFsteel') ## moduls:
abaquspy.mats.ifsteel.iso(myMat)

# ## Create Shell Section!
thickness = 1e-3 ## 1 mm thickness
myShellSection=myModel.HomogeneousShellSection(
    name='SpecimenSection',
    preIntegrate=OFF, material=myMat.name, thicknessType=UNIFORM, thickness=thickness,
    thicknessField='', idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT,
    thicknessModulus=None, temperature=GRADIENT, useDensity=OFF,
    integrationRule=SIMPSON, numIntPts=5)

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
totalStrain = 0.1
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
# myModel.historyOutputRequests['H-Output-1'].setValues(variables=('E11',),region=myAssembly.sets['MidSpan'])


# ## Apply BC
c0=datOri.pointOn
vs=myInstance.vertices.findAt((c0,))
## Encastre
myModel.EncastreBC(name='EncastreOri',createStepName='Tension',region=regionToolset.Region(vertices=vs))
## Symmetric constraints
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
elemType2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
#f = myPart.faces
faces=myPart.faces[:]# = f.getSequenceFromMask(mask=('[#f ]', ), )
pickedRegions =(faces, )
myPart.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
myPart.seedPart(size=0.005, minSizeFactor=0.1) ## coarse meshing
#myPart.seedPart(size=0.001, minSizeFactor=0.1) ## finer meshing
myPart.generateMesh()

## additional procedure that requires 'meshed' information.
myAssembly.regenerate()
myAssembly.Set(elements=myInstance.elements[:],name='ORIGIN')
myAssembly.regenerate()

myModel.HistoryOutputRequest(name='StressStrain',
                             createStepName='Tension',variables=('S11','E11','E22','PE11','PE22'),
                             region=myAssembly.sets['ORIGIN'],sectionPoints=DEFAULT,rebar=EXCLUDE)
myAssembly.regenerate()

## Create Job
mdb.Job(name='OneElement',model=myModel.name,description='PythonScriptedOneElement')
mdb.saveAs(myModel.name)
myJob = mdb.jobs['OneElement']
## Multicore options
myJob.setValues(numCpus=4,numDomains=4)


## Flag to use a User Material subroutine
iumat=False
umatFN='/home/younguj/repo/abaqusPy/umats/el/iso.f'
if iumat:
    myJob.setValues(userSubroutine=umatFN)

if False:
    ## submit the job
    myJob.submit(consistencyChecking=OFF)
    ## myJob wait until completion?
    myJob.waitForCompletion()
    ## execute pp file?
    # execfile('onePP.py')

    #import onePP
    #onePP.main(session,'strstr_%.2f.txt'%Theta)
