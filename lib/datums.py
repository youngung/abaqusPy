import numpy as np
from abaqusConstants import *


def returnPointDatumCoords(myPart,coords):
    """
    Define a point datum using the given coordinates

    Arguments
    --------
    coords
    """

    if type(coords)==tuple: pass
    elif type(coords)==np.ndarray or type(coords)==list:
        coords=tuple(coords)
    else:
        raise SyntaxError,'Unexpected type of coords given'

    ## Below returns a feature
    f=myPart.DatumPointByCoordinate(coords=coords)
    ## Convert the feature f to the datum object using its id
    return myPart.datums[f.id]

def returnCsymPlanarAngle(myPart,datOri,radius,angle):
    """
    This function defines two point datums
    then create cooridate system datum that is
    in the plane of XY of the global axes and is
    in c.c.w. rotation from X inasmuchas <angle> in radian.


    Arguments
    ---------
    myPart  : Abaqus Part object
    datOri  : datum that gives original orientation]
    radius  : Radius of the coordinate axes  [in meters]
    angle   : angle inclined ccw from global x direction   [in Radia]
    """
    import numpy as np
  ## define the coordiantes for two point datums.

    cOri = np.array(datOri.pointOn)

    xvector = np.array([np.cos(angle), np.sin(angle),0])*radius
    yvector = np.array([np.cos(angle+np.pi/2.), np.sin(angle+np.pi/2.),0])*radius

    x = cOri + xvector
    y = cOri + yvector

    datX = returnPointDatumCoords(myPart,tuple(x))
    datY = returnPointDatumCoords(myPart,tuple(y))

    cSysMat = myPart.DatumCsysByThreePoints(
        CARTESIAN,origin=datOri,
        point1=datX,point2=datY,name='matCsys')

    return cSysMat
