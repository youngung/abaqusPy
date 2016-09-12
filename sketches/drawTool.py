def draw_arc(mySketch,xyCoords,radius,iref,i0,i1,opt='up',):
    """
    arguments
    ---------
    mySketch: Abaqus sketch object
    xyCoords: collection of xy coordinates
    radius
    iref: index of coordinates to which circle center is referenced.
    i0  : arc start index.
    i1  : arc end index
    opt : 'up' or 'down' depending on the reletive location of circle center to specimen.
    """
    xyArc=xyCoords[iref]

    if opt=='up':    centerArc = (xyArc[0],xyArc[1]+radius)
    elif opt=='down':centerArc = (xyArc[0],xyArc[1]-radius)
    else: raise IOError, 'Unexpected opt given'

    p1=(xyCoords[i0][0],xyCoords[i0][1])
    p2=(xyCoords[i1][0],xyCoords[i1][1])
    mySketch.ArcByCenterEnds(center=centerArc, point1=p1, point2=p2)

def draw_line(mySketch,xyCoords,i0,i1):
    """
    Arguments
    ---------
    xyCoords
    i0
    i1
    """
    n=i1-i0
    for i in xrange(n):
        I0=i0+i
        I1=I0+1
        p1=(xyCoords[I0][0],xyCoords[I0][1])
        p2=(xyCoords[I1][0],xyCoords[I1][1])
        mySketch.Line(point1=p1,point2=p2)
    pass
