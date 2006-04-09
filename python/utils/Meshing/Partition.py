#!/usr/bin/env python

import Chaco
from Chaco import *

import TriangleReader
from TriangleReader import *

import meshUtils
from meshUtils import *

offset = 0

r = TriangleReader('diskWithHole-verycoarse', offset)

mesh = r.getMesh()


filename = 'diskWithHole-verycoarse'

c = Chaco(filename, offset)

np = 4

elemAssignments = c.partition(mesh, np)


(nodeAssignments, nodeOwners, nodesPerProc) = mesh.getNodeAssignments(np, elemAssignments)

elemsPerProc = mesh.getElemsPerProc(np, elemAssignments)



for procID in range(np) :
    (offProcNodes, offProcElems) = mesh.getOffProcData(procID, elemAssignments,
                                                       nodeAssignments)
    nodeGIDToLIDMap = writeNodes(filename, mesh, procID, np, nodeAssignments,
                                 nodesPerProc[procID], offProcNodes)
    writeElems(filename, mesh, procID, np, elemAssignments,
               elemsPerProc[procID], offProcElems, nodeGIDToLIDMap)

    (nodeGID, nodeOwners) = mesh.getNodeParInfo(procID, np, nodeAssignments, offProcNodes)
    (elemGID, elemOwners) = mesh.getElemParInfo(procID, np, elemAssignments, offProcElems)
    writeParFile(filename, procID, np, nodeGID, nodeOwners, elemGID, elemOwners)
    
