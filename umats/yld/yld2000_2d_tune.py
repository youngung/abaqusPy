import numpy as np

def virtual_test(alphas,exponent):
    """
    Conduct virtual uniaxial tension tests (and balanced biaxial)
    return the 'difference'

    Arguments
    ---------
    alphas
    exponent
    """
    from yld2000 import yld2000_2d, inplane_rot, voigt1,voigt2,\
        voigt4,reduce_6to3, reduce_3to6

    ## alphas and exponent
    yldc = np.zeros(9)
    yldc[:8] = alphas
    yldc[8]  = exponent

    ## uniaxial directions to probe
    nth = 3
    ths = np.linspace(0,np.pi/2.,nth)
    ys = np.zeros(nth+1)
    rs = np.zeros(nth+1)

    ## 6D stress state
    s6lab = np.zeros(6)
    s6lab[:] = 0
    s6lab[0] = 1.
    s33lab = voigt2(s6lab)

    for i in range(len(ths)):
        ## stress in material coordinates
        s33mat = inplane_rot(ths[i],s33lab)
        ## convert to 6D convenction
        s6mat = voigt1(s33mat)
        s3mat = reduce_6to3(s6mat)
        phi,dphi,d2phi = yld2000_2d(s3mat,yldc)
        dphi = reduce_3to6(dphi)
        ## convert dphi to dphi33 in mat
        dphi33m = voigt4(dphi)
        ## convert dphi in mat to dphi in lab
        dphi33l = inplane_rot(-ths[i],dphi33m)
        ## r-value
        rv = -dphi33l[1,1] / (dphi33l[0,0]+dphi33l[1,1])
        rs[i]=rv
        ## yield stress
        ys[i]=phi

    ## conduct balanced biaxial
    s6mat      = np.zeros(6)
    s6mat[0:2] = 1.0
    s3mat = reduce_6to3(s6mat)
    ys_bb,dphi,d2phi = yld2000_2d(s3mat,yldc)
    rv_bb = dphi[1]/dphi[0]

    ys[-1] = ys_bb
    rs[-1] = rv_bb

    return ys,rs

def diff_constant_exponent(alphas,exponent,rs_exp,ys_exp):
    """
    Arguments
    ---------
    alphas
    exponent
    rs_exp
    ys_exp
    """
    ys_mod, rs_mod = virtual_test(alphas,exponent)
    f = (ys_exp/ys_exp[0] - ys_mod)**2 + (rs_mod - rs_exp)**2
    return f.sum()

def diff_wrapper_constant(exponent,r0,r45,r90,rb,y0,y45,y90,yb):
    rs_exp = np.array([r0,r45,r90,rb])
    ys_exp = np.array([y0,y45,y90,yb])

    def func(alphas):
        return diff_constant_exponent(alphas,exponent,rs_exp,ys_exp)
    return func

def diff(yldc,rs_exp,ys_exp):
    alphas=yldc[:-1]
    exponent=yldc[-1]
    ys_mod, rs_mod = virtual_test(alphas,exponent)
    f = (ys_exp/ys_exp[0] - ys_mod)**2 + (rs_mod - rs_exp)**2
    return f.sum()

def diff_wrapper(r0,r45,r90,rb,y0,y45,y90,yb):
    rs_exp = np.array([r0,r45,r90,rb])
    ys_exp = np.array([y0,y45,y90,yb])
    def func(yldc):
        return diff(yldc,rs_exp,ys_exp)
    return func

def tune1(*yldc):
    """
    VPSC tunning using yldc that includes exponent
    The last element in given yldc is regarded as a
    prescribed exponent to yield surface
    """
    import scipy.optimize
    func = diff_wrapper(*yldc)

    x0=np.ones(9)
    x0[-1]=2

    xport = scipy.optimize.fmin(
        func=func,x0=x0)
    return xport

def tune2(exponent=8,*yldc):
    import scipy.optimize
    func = diff_wrapper_constant(exponent,*yldc)
    x0=np.ones(8)
    xport = scipy.optimize.fmin(
        func=func,x0=x0)
    return xport
