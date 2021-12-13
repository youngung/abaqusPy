"""
Python script to analyze the resulting data from uni.py

Youngung Jeong
youngung.jeong@gmail.com
"""

from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *

fn='/home/younguj/abaqus/exercise/uten/TensileE8.odb'
odb = session.openOdb(name=fn)

## Analyze the model
print('Instances in the assembly')
for instanceName in list(odb.rootAssembly.instances.keys()):
    print(instanceName)

# - material
allMaterials = odb.materials
for materialName in list(allMaterials.keys()):
    print('Material Name : ',materialName)
#To print isotropic elastic material properties in a material object:

for material in list(allMaterials.values()):
    if hasattr(material,'elastic'):
        elastic = material.elastic
        if elastic.type == ISOTROPIC:
            print('isotropic elastic behavior, type = %s' \
            % elastic.moduli)
        title1 = '%9s %6s %9s'%('Young\'s modulus','','Poisson ratio')
        title2 = ''
        if elastic.temperatureDependency == ON:
            title2 = 'Temperature  '
        dep = elastic.dependencies
        title3 = ''
        for x in range(dep):
            title3 += ' field # %d' % x
        print('%s %s %s' % (title1,title2,title3))
        for dataline in elastic.table:
            y=dataline[0]
            nu=dataline[1]
            print('%9.1f %6s %9.1f'%(y*1e-9,'[GPa]',nu))

session.Viewport(name='Viewport: 2')            
            
session.viewports['Viewport: 2'].setValues(displayedObject=odb)
odbName=session.viewports['Viewport: 2'].odbDisplay.name

print('%5s %20s %20s'%('Id','Step Name','Description'))
stepContainer = list(odb.steps.keys())
for key in list(odb.steps.keys()):
    step=odb.steps[key]
    print('%5i %20s %20s'%(step.number,step.name, step.description))

stepObj=odb.steps[stepContainer[-1]] ## deal with the last step container only.

## count how many 'frames' exist?
lastFrame=stepObj.frames[-1]
nFrame=lastFrame.frameId


# total number
TotalNumber=lastFrame.fieldOutputs['E'].values.__len__()
lastFrame.fieldOutputs['E'].values[0].data[0]## strain Exx?



## data from RIGHTEND set
setRightEnd=odb.rootAssembly.instances['MYSPECIMEN'].nodeSets['RIGHTEND']

## integration point
#StrainField = lastFrame.fieldOutputs['E'].getSubset(region=setRightEnd,position=INTEGRATION_POINT,elementType='S4R')
## nodal
StrainField = lastFrame.fieldOutputs['E'].getSubset(region=setRightEnd,position=NODAL,elementType='S4R')

j=0
for v in StrainField.values:
    v.data ## data array
    j=j+1
    for i in range(len(v.data)):
        print('component %i: %.3f'%(i+1,v.data[i]))
        
print('j:',j)









