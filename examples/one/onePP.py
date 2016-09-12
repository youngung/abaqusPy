"""
Conduct analysis on OneElement.odb file
"""
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import os

import numpy as np

def main(session,fnout='strstr.txt'):
    """

    Arguments
    ---------
    session  : abaqus session object
    fnout    : output filename
    """
    fnOdb=os.path.join(os.getcwd(),'OneElement.odb')
    print 'fnOdb:',fnOdb
    odb = session.openOdb(name=fnOdb)

    ## Analyze the model
    print 'Instances in the assembly'
    for instanceName in odb.rootAssembly.instances.keys():
        print instanceName

    # - material
    allMaterials = odb.materials
    for materialName in allMaterials.keys():
        print 'Material Name : ',materialName
    #To print isotropic elastic material properties in a material object:

    for material in allMaterials.values():
        if hasattr(material,'elastic'):
            elastic = material.elastic
            if elastic.type == ISOTROPIC:
                print 'isotropic elastic behavior, type = %s' \
                % elastic.moduli
            title1 = '%9s %6s %9s'%('Young\'s modulus','','Poisson ratio')
            title2 = ''
            if elastic.temperatureDependency == ON:
                title2 = 'Temperature  '
            dep = elastic.dependencies
            title3 = ''
            for x in range(dep):
                title3 += ' field # %d' % x
            print '%s %s %s' % (title1,title2,title3)
            for dataline in elastic.table:
                y=dataline[0]
                nu=dataline[1]
                print '%9.1f %6s %9.1f'%(y*1e-9,'[GPa]',nu)

    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    odbName=session.viewports['Viewport: 1'].odbDisplay.name

    print '%5s %20s %20s'%('Id','Step Name','Description')
    stepContainer = odb.steps.keys()
    for key in odb.steps.keys():
        step=odb.steps[key]
        print '%5i %20s %20s'%(step.number,step.name, step.description)

    stepObj=odb.steps[stepContainer[-1]] ## deal with the last step container only.

    ## count how many 'frames' exist?
    lastFrame=stepObj.frames[-1]
    nFrame=lastFrame.frameId

    ## extract E11/S11
    #odb = session.odbs['/home/younguj/abaqus/exercise/uten/one/OneElement.odb']

    ## total strain
    e11=odb.steps['TensionContinue'].historyRegions[
        'Element MYSPECIMEN.1 Int Point 1 Section Point 1'].historyOutputs['E11'].data
    pe11=odb.steps['TensionContinue'].historyRegions[
        'Element MYSPECIMEN.1 Int Point 1 Section Point 1'].historyOutputs['PE11'].data
    pe22=odb.steps['TensionContinue'].historyRegions[
        'Element MYSPECIMEN.1 Int Point 1 Section Point 1'].historyOutputs['PE22'].data
    s11=odb.steps['TensionContinue'].historyRegions[
        'Element MYSPECIMEN.1 Int Point 1 Section Point 1'].historyOutputs['S11'].data
    e11=np.array(e11);    s11=np.array(s11)
    pe11=np.array(pe11); pe22=np.array(pe22)
    ## Discard the time stamps.
    time=e11[:,0]
    s11=s11[:,1]
    e11=e11[:,1];
    pe11=pe11[:,1];
    pe22=pe22[:,1]; 

    FlowCurve=np.array([e11,s11,pe11,pe22,time]).T
    fnFlowCurve=os.path.join(fnout)
    np.savetxt(fnFlowCurve,FlowCurve)
