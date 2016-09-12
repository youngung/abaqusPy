"""
Interstital-free steel
"""
import os
os.sys.path.append('/home/younguj/anaconda2/lib/python2.7/site-packages/')
import mk.materials.func_hard as func_hard
import numpy as np

def hard(iopt=0,emx=0.1,nstp=1000):
    if iopt==0:
        ## use voce fit?
        ## Voce BB fit [MPa]
        params='403.23447661  273.35755903   11.95412158  191.28660496'
        params=map(float,params.split())

    elif iopt==1:
        ## Voce UR fit [MPa]
        params='255.51670864  129.97608057   28.98133641  454.07524122'
        params=map(float,params.split())
    else:
        raise SyntaxError, 'Unexpected iopt given'

    eps=np.linspace(0,emx,nstp)
    sig=func_hard.func_voce(eps,*params)
    sig=sig*1e6 ## to Pascal

    # make the data to tuple
    dat=np.array([sig,eps]).T
    return tuple(map(tuple,dat))

def elastic_iso(myMat):
    """
    Isotropic elasticity

    Argument
    --------
    """
    young = 200e9 ## young's modulus
    nu    = 0.3   ## poisson ratio
    myMat.Elastic(table=((young,nu),))

def plastic_iso(myMat):
    """
    Argument
    --------
    myMat
    """
    hardDatTable = hard(iopt=0,emx=0.5,nstp=100)
    myMat.Plastic(table=hardDatTable[::])


## isotropic case...
def iso(myMat):
    """
    """
    elastic_iso(myMat) ## isotropic elastic
    plastic_iso(myMat) ## isotorpic plastic
