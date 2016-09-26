"""
Conduct analysis on OneElement.odb file
"""
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import os

import numpy as np

def main(fnout='strstr.txt',fnOdb=None,iumat=False):
    """

    Arguments
    ---------
    fnout    : output filename
    iumat    : False (flag to determine if user material is used.
    """
    ## fnOdb=os.path.join(os.getcwd(),'OneElement.odb')
    print 'fnOdb:',fnOdb
    odb = openOdb(fnOdb)

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
    ## total strain
    if odb.rootAssembly.instances.__len__()>1:
        print 'warning: expected single instance but found multiple'
    if odb.sections.__len__()>1:
        print 'warning: expected single section but found multiple'

    myInstancName = odb.rootAssembly.instances.keys()[0]
    mySectionName = odb.sections.keys()[0]

    hRegion='Element MYSPECIMEN.1 Int Point 1 Section Point 1'



    ## In case intrinsic material feature in abaqus is used.
    if not(iumat):
        ## Commonly available variables
        e11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['E11'].data
        s11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['S11'].data
        ee11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['EE11'].data
        ee22=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['EE22'].data
        pe11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['PE11'].data
        pe22=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['PE22'].data
    elif (iumat):
        # ee11=odb.steps['TensionContinue'].historyRegions[hRegion].\
        #     historyOutputs['SDV1'].data
        # ee22=odb.steps['TensionContinue'].historyRegions[hRegion].\
        #     historyOutputs['SDV2'].data
        # pe11=odb.steps['TensionContinue'].historyRegions[hRegion].\
        #     historyOutputs['SDV4'].data
        # pe22=odb.steps['TensionContinue'].historyRegions[hRegion].\
        #     historyOutputs['SDV5'].data
        ee11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['UVARM1'].data
        ee22=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['UVARM2'].data
        pe11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['UVARM4'].data
        pe22=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['UVARM5'].data
        s11=odb.steps['TensionContinue'].historyRegions[hRegion].\
            historyOutputs['UVARM7'].data
        ee11=np.array(ee11)
        ee22=np.array(ee22)
        pe11=np.array(pe11)
        pe22=np.array(pe22)
        e11=ee11+pe11
        e22=ee22+pe22
    else:
        print 'Unexpected iumat given:',iumat
    ## In case intrinsic material feature in abaqus is used.


    e11=np.array(e11);    s11=np.array(s11)
    pe11=np.array(pe11); pe22=np.array(pe22)
    ee11=np.array(ee11); ee22=np.array(ee22)
    ## Discard the time stamps.
    time=e11[:,0]
    s11=s11[:,1]
    e11=e11[:,1]
    pe11=pe11[:,1]
    pe22=pe22[:,1]
    ee11=ee11[:,1]
    ee22=ee22[:,1]

    FlowCurve=np.array([e11,s11,pe11,pe22,time,ee11,ee22]).T
    fnFlowCurve=os.path.join(os.getcwd(),fnout)
    np.savetxt(fnFlowCurve,FlowCurve,fmt='%13.5e')

    print 'fnFlowCurve file %s has been saved using onePP.py'%fnFlowCurve

    return odb

print '\n\n\n'

import glob
fns=glob.glob('*.odb')
#(fnout='strstr.txt',fnOdb=None
for i in xrange(len(fns)):
    print 'Concerned file name:',fns[i]
    odb=main(fnout=fns[i].split('.odb')[0]+'.txt',fnOdb=fns[i],iumat=True)
