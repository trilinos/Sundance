#! /usr/bin/env python


import PySundance
import math
from PySundance import *


#############################################################################
#
# Read a Triangle file, then write it to VTK
#
#############################################################################

def main():

    # Read the mesh from a Triangle file
    filename = 'diskWithHole-med'
    mesher  = TriangleMeshReader(filename)
    mesh = mesher.getMesh()
    check = mesh.checkConsistency('meshCheck')
    print 'check=', check
    if check==0 :
        print 'INCONSISTENT MESH'
    
    vecType = EpetraVectorType()
    basis = Lagrange(0)
    discSpace = DiscreteSpace(mesh, basis, vecType)

    p = float(getRank()) / float(getNProc())
    f = DiscreteFunction(discSpace, p)

    writer = VTKWriter(filename)
    writer.addMesh(mesh)
    writer.addField('procID', f)
    writer.write()

    # all done!

    

if __name__ == "__main__":
    main()

    







    
    
