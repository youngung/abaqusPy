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
from driverUtils import executeOnCaeStartup
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

def fetch_input(name='springback_std_both'):
    cmd='abaqus fetch job=%s'%name
    os.system(cmd)


## reference model using fetched input file from ABAQUS
modelName_reference='springback_std_both'
fetch_input(modelName_reference)
## executeOnCaeStartup()
model_reference = mdb.ModelFromInputFile(name=modelName_reference,inputFileName='%s.inp'%modelName_reference)


## New model file by copying the reference file to
modelName='%s_umat'%modelName_reference
myModel = mdb.Model(name=modelName,
                   objectToCopy=model_reference)


myModel.fieldOutputRequests['F-Output-1'].setValues(frequency=1)
myModel.fieldOutputRequests['F-Output-2'].setValues(frequency=1)
myModel.fieldOutputRequests['F-Output-2'].setValues(variables=('S', 'E','UVARM'))

myPart=myModel.parts['PART-1']

## orientation
datFeature = myPart.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
        0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
region = myPart.sets['BLANK']
orientation = myPart.datums[3]
myPart.MaterialOrientation(
    region=region, orientationType=SYSTEM, axis=AXIS_3, localCsys=orientation,
    fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0,
    additionalRotationField='')


myMat = myModel.Material('myUMAT')
myAssembly = myModel.rootAssembly
myMat.UserMaterial(mechanicalConstants=(0.0,))
myMat.Depvar(n=20)
myMat.UserOutputVariables(n=20) ## UVARM

gpa=1e9
mpa=1e6

myAssembly.regenerate()
myModel.sections['Section-1-BLANK'].setValues(
    preIntegrate=OFF, material='myUMAT', thicknessType=UNIFORM,
    thickness=0.00078, thicknessField='', idealization=NO_IDEALIZATION,
    integrationRule=SIMPSON, numIntPts=5)
myModel.sections['Section-1-BLANK'].TransverseShearShell(
    k11=2.*mpa,k22=2.*mpa,k12=1.2*mpa)
myAssembly.regenerate()

jobName='springback_std_both_umat'
mdb.Job(name=jobName,model=myModel.name,
        description='PythonScripted_%s'%jobName)
jobName_ref='springback_std_both'
mdb.Job(name=jobName_ref,model=myModel.name,
        description='PythonScripted_%s'%jobName)

mdb.saveAs(myModel.name)
myJob = mdb.jobs[jobName]
mdb.saveAs(model_reference.name)
myJob_ref = mdb.jobs[jobName_ref]

## Flag to use a User Material subroutine
print 'User material has been specified.'
umatFN='/home/younguj/repo/abaqusPy/umats/epl/epl.f'
myJob.setValues(userSubroutine=umatFN)
myJob.submit(consistencyChecking=OFF)
myJob.waitForCompletion()
