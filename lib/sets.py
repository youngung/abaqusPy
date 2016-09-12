"""
Function to define 'sets'
"""
def setNodeCoord(myPart,dat=None,name='Nodal-Set-1',offset=1e-4):
    """
    Define a set attribute within myInstance using
    the given datum <dat>. Select a <single> node
    <datum> is assumed to be a point datum.

    Arguments
    ---------
    myPart - Abaqus part object to which the set will be assinged
    dat    - this datum shuould be a point datum.
    name   - Name for the new set.
    offset - coordinate offset
    """
    import numpy as np
    if type(dat).__name__!='DatumPoint':
        raise SyntaxError, 'Datum should be DatumPoint type'

    coordinate = np.array(dat.pointOn)
    x,y,z=tuple(coordinate)
    xmin,ymin,zmin=tuple(coordinate-offset)
    xmax,ymax,zmax=tuple(coordinate+offset)

    ## assign node to "part". Regenerating assembly will inherit the node.
    selectedNodes=myPart.nodes.getByBoundingBox(xmin,ymin,zmin,xmax,ymax,zmax)
    myPart.Set(nodes=(selectedNodes),name=name)
