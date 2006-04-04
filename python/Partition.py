#!/usr/bin/env python

import Chaco
from Chaco import *

import TriangleReader
from TriangleReader import *

import meshUtils
from meshUtils import *

offset = 0

r = TriangleReader('halfCyl', offset)

mesh = r.getMesh()

c = Chaco('joe', offset)

np = 128

elemAssignments = c.partition(mesh, np)


(nodeAssignments, nodeOwners, nodesPerProc) = mesh.getNodeAssignments(np, elemAssignments)

elemsPerProc = mesh.getElemsPerProc(np, elemAssignments)

distributeNodes('joe', mesh, np, nodeAssignments, nodeOwners, nodesPerProc)

distributeElems('joe', mesh, np, elemAssignments, elemsPerProc)
