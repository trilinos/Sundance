############################################################################
#
# Python wrapper for Metis partitioning 
#
############################################################################

from sets import Set
from Mesh import Mesh
import posix
import string

class Metis :

    # create a partitioner object
    def __init__(self, filename) :
        self.filename_ = filename


    # driver routine for the partitioning
    def partition(self, mesh, nProc) :

        mesh.writeGraph(self.filename_)
        self.runMetis(nProc)
        return self.readElemAssignments(self.filename_)

    # Set up metis input file and run metis

    def runMetis(self, nProcs) :
        posix.system('kmetis %s.graph %d' % (self.filename_, nProcs) )
        posix.system('cp %s.graph.part.%d %s.assign' % (self.filename_, nProcs, self.filename_) )

    def readElemAssignments(self, filename) :
        assignments = []
        f = file(self.filename_ + '.assign')

        # read lines in the assignment file
        while 1:
            line = f.readline()
            if not line : break
            if line[0]=='#': continue
            assignments.append(int(line))

        return assignments


    def writePartitionFile(self, filename, assignments, nProcs) :
        partFile = file(filename + '.part', 'w')
        partFile.write('%d %d\n' % (len(assignments), nProcs))
        for i in range(len(assignments)) :
            partFile.write('%d %d\n' % (i, assignments[i]))
