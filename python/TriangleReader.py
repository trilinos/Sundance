#@HEADER@
#@HEADER@


############################################################################
#
# Read Triangle file and create Python mesh object
#
############################################################################

from Mesh import Mesh
from sets import Set

class TriangleReader :

    def __init__(self, filename, indexOffset) :
        self.filename_ = filename
        self.indexOffset_ = indexOffset

    def getMesh(self) :

        mesh = Mesh()

        self.readPoints(mesh)

        self.readElems(mesh)

        return mesh

    # Read the .node file
    def readPoints(self, mesh) :

        f = file(self.filename_ + '.node')

        # read header, which looks like [nNodes, dim, nAttr, nBdry] 
        while 1 :
            line = f.readline()
            if line[0]=='#': continue
            headerline = line
            header = line.split()
            nNodes = int(header[0])
            d = int(header[1])
            nAttr = int(header[2])
            nMark = int(header[3])
            break

        mesh.setDimension(d)

        # read remaining lines, adding the points
        # each line looks like [gid, x_1 ... x_d , attrs, bdry]
        while 1:
            line = f.readline()
            if not line : break
            if line[0]=='#': continue
            data = map(float, line.split()[1:d+1])
            mesh.addPoint(data)


    # Read the .ele file
    def readElems(self, mesh) :
        
        f = file(self.filename_ + '.ele')

        # read header 
        while 1 :
            line = f.readline()
            if line[0]=='#': continue
            header = line.split()
            nElems = int(header[0])
            d = int(header[1])-1
            break

         # read lines, building elements and the element-to-node map
        while 1:
             line = f.readline()
             if not line : break
             if line[0]=='#': continue
             toks = line.split()
             ele = int(toks[0])
             verts = Set()
             for i in range(d+1) :
                 node = int(toks[i+1])-self.indexOffset_
                 verts.add(node)
             mesh.addElem(verts)
        
    
    
