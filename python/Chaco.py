#@HEADER@
#@HEADER@

############################################################################
#
# Python wrapper for Chaco partitioning 
#
############################################################################

from sets import Set
from Mesh import Mesh
import posix
import string

class Chaco :

    # create a partitioner object
    def __init__(self, filename, indexOffset) :
        self.filename_ = filename
        self.indexOffset_ = indexOffset


    # driver routine for the partitioning
    def partition(self, mesh, nProc) :

        mesh.writeGraph(self.filename_)
        self.runChaco(nProc)
        return self.readElemAssignments(self.filename_)

    # Set up chaco input file and run chaco

    def runChaco(self, nProcs) :
        paramsFile = file('User_Params', 'w')
        paramsFile.write('OUTPUT_ASSIGN=true\n')
        paramsFile.write('PROMPT=false\n')
        paramsFile.write('ARCHITECTURE=1\n')
        paramsFile.write('REFINE_PARTITION=4\n')
        paramsFile.write('REFINE_MAP=true\n')
        paramsFile.write('KL_BAD_MOVES=20\n')
        paramsFile.write('KL_NTRIES_BAD=10\n')
        paramsFile.write('KL_IMBALANCE=0.02\n')
        paramsFile.write('INTERNAL_VERTICES=true\n')
        paramsFile.write('MATCH_TYPE=4\n')
        paramsFile.write('HEAVY_MATCH=true\n')
        paramsFile.write('TERM_PROP=true\n')
        paramsFile.write('COARSE_NLEVEL_KL=1\n')
        paramsFile.write('COARSEN_RATIO_MIN=0.7\n')
        paramsFile.write('CUT_TO_HOP_COST=1.0\n')
        paramsFile.write('RANDOM_SEED=12345\n')
        paramsFile.flush()
        inputFile = file('chacoInput', 'w')
        inputFile.write('%s.graph\n' % self.filename_)
        inputFile.write('%s.assign\n' % self.filename_)
        inputFile.write('1\n')
        inputFile.write('100\n')
        inputFile.write('%d\n' % nProcs)
        inputFile.write('1\n')
        inputFile.write('n\n')
        inputFile.flush()
        inputFile.close()
        posix.system('chaco < chacoInput')

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
            partFile.write('%d %d\n' % (i+self.indexOffset_, assignments[i]+self.indexOffset_))
