"""
modelAExample.py

A simple example: Creating a part.
"""

import os

from abaqus import *
from abaqusConstants import *
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

import sketch
import part
import numpy as np


###### ----------Generating specimen dimensions.

def drawE8_Sketch(
    pl=57,   ## Parallel length
    gw=12.5, ## Gauge width
    tw=20,   ## Grip width
    tl=50.0, ## Grip length
    rd=12.5, ## Radius
    ):
    """
    Draw E8 uniaxial tension specimen
    Units are in [mm]
    
    
    Arguments
    ---------
    pl : parallel length  (57   mm)
    gw : gauge width      (12.5 mm)
    tw : grip width       (20   mm)
    tl : grip length      (50.0 mm)
    rd : radius           (12.5 mm)
    """

    A  = tw/2. - gw/2.
    th = -90.-np.arccos((rd-A)/rd)*180./np.pi
    x  = rd * np.sin(th*np.pi/180.)

    ## Round...
    th0=-90
    th_delta=np.arccos((rd-A)/rd)*180./np.pi
    th1=th0+th_delta
    # ths=np.linspace(th0*np.pi/180.,th1*np.pi/180.)
    ths=[th0*np.pi/180.,th1*np.pi/180.]
    xs = rd*np.cos(ths)
    ys = rd*np.sin(ths)

    ## translate xs,ys
    xs = xs + (x-xs[-1])
    ys = ys + (-A+tw-ys[-1])
    xyRound=[xs.tolist(),ys.tolist()]
    
    ## parallel
    x0,y0=xs[0],ys[0]
    xyParallel = [[x0-0.5*pl,x0],[y0,y0]]
    
    ## Right grip
    XS=[x+tl,x+tl,x][::-1]
    YS=[-A+0.5*tw,-A+tw,-A+tw][::-1]
    xyRG=[XS,YS]
    
    x=xyParallel[0]+xyRound[0]+xyRG[0]
    y=xyParallel[1]+xyRound[1]+xyRG[1]
    
    xyCoords=np.array([x,y])
    
    # print xyCoords.shape
    
    ## translate the coordinate so that the center of gravity is (0,0)
    xyCoords[0]=xyCoords[0]-xyCoords[0][0]
    xyCoords[1]=xyCoords[1]-xyCoords[1][-1]
    # plot(xyCoords[0],xyCoords[1],'-')
    
    ## Apply 2-fold symmetry.
    sym0 =[[ 1,0],[0, 1]] ## Identical
    sym1 =[[-1,0],[0, 1]] ## Reflect y axis
    sym2 =[[ 1,0],[0,-1]] ## Reflect x axis
    sym3 =[[-1,0],[0,-1]] ## Reflect origin
    
    sym = np.array([sym0,sym2,sym3,sym1])
    # plot(xyCoords[0,0],xyCoords[1,0],'x')

    xyTot=[[],[]]
    for i in xrange(len(sym)):
        symOp = sym[i][:,:]# (2,2)
        temp = np.tensordot(symOp,xyCoords,axes=[1,0])
        if i==1 or i==3:
            temp[0]=temp[0][::-1]
            temp[1]=temp[1][::-1]
        elif i==0 or i==2:
            temp=temp[::]

        for j in xrange(len(temp[0])):
            xyTot[0].append(temp[0][j])
            xyTot[1].append(temp[1][j])

    xyTot=np.array(xyTot)


    x0=min(xyTot[0])
    y0=min(xyTot[1])+tw/2.
    
    xyTot[0] = xyTot[0] - x0
    xyTot[1] = xyTot[1] - y0
    
    return xyTot

# ---
def tensileBar(pl=57,   ## Parallel length
               gw=12.5, ## Gauge width
               tw=20,   ## Grip width
               tl=50.0, ## Grip length
               rd=12.5, ## Radius
               ):                   
    """
    Close end.
    """
    xy=drawE8_Sketch(pl,gw,tw,tl,rd)
    xy[1,:] = xy[1,:]+tw/2.
    xyt=xy.T
    x,y=xy
    XYS=[]
    for i in xrange(len(xyt)-1):
        x0,y0 = xyt[i]
        x1,y1 = xyt[i+1]
        if (x1-x0)==0 and (y1-y0)==0:
            pass
        else:
            XYS.append([x1,y1])
            
    XYS.append(XYS[0])
    xyt=np.array(XYS).T
    return xyt
## --

## dimension in meter
pl=57.    * 1e-3
gw=12.5   * 1e-3
tw=20.    * 1e-3
tl=50.0   * 1e-3
rd=12.5   * 1e-3 ## radius
xyCoords = tensileBar(pl=pl,gw=gw,tw=tw,tl=tl,rd=rd).T

myModel = mdb.Model(name='UniaxialTension')
mySketch = myModel.ConstrainedSketch(
    name='E8',sheetSize=1.0)
    
#s.ArcByCenterEnds(center=?, point1=?,point2?)

for i in xrange(len(xyCoords)):
    mySketch.Spot(point=xyCoords[i])

#for i in xrange(len(xyCoords)-1):
#    mySketch.Line(point1=(xyCoords[i,0],xyCoords[i,1]),
#                  point2=(xyCoords[i+1,0],xyCoords[i+1,1]))
                  
g,v,d,c=mySketch.geometry, mySketch.vertices, mySketch.dimensions, mySketch.constraints

def draw_arc(iref,i0,i1,opt='up',):
    """
    arguments
    ---------
    iref: index of coordinates to which circel center is referenced.
    i0  : arc start index.
    i1  : arc end index
    opt : 'up' or 'down' depending on the reletive location of circle center to specimen.
    """
    xyArc=xyCoords[iref]
    
    if opt=='up':    centerArc = (xyArc[0],xyArc[1]+gw)
    elif opt=='down':centerArc = (xyArc[0],xyArc[1]-gw)
    else: raise IOError, 'Unexpected opt given'

    p1=(xyCoords[i0][0],xyCoords[i0][1])
    p2=(xyCoords[i1][0],xyCoords[i1][1])
    mySketch.ArcByCenterEnds(center=centerArc, point1=p1, point2=p2)
    
def draw_line(i0,i1):
    n=i1-i0
    for i in xrange(n):
        I0=i0+i
        I1=I0+1
        p1=(xyCoords[I0][0],xyCoords[I0][1])
        p2=(xyCoords[I1][0],xyCoords[I1][1])
        
        mySketch.Line(point1=p1,point2=p2)

draw_arc(0,0,1,'up')
draw_arc(8,7,8,'down')
draw_arc(10,10,11,'down')
draw_arc(18,17,18,'up')

draw_line(1,7)
draw_line(8,10)
draw_line(11,17)
draw_line(18,20)


###### ----------Generating specimen dimensions.



myPart = myModel.Part(name='E8', dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)
