/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceMap.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceDOFMapBase.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;




DOFMapBase::DOFMapBase(const Mesh& mesh)
  : localProcID_(mesh.comm().getRank()),
    mesh_(mesh),
    lowestLocalDOF_(),
    numLocalDOFs_(),
    numDOFs_(),
    ghostIndices_(rcp(new Array<int>()))
{}

void DOFMapBase::verifySubregionIndex(int subregionIndex) const
{
  TEST_FOR_EXCEPTION(subregionIndex < 0 
                     || subregionIndex >= numHomogeneousSubregions(),
                     RuntimeError,
                     "Subregion index=" << subregionIndex
                     << " out of range [0, " << numHomogeneousSubregions()
                     << ")");
}

void DOFMapBase::verifyBasisChunkIndex(int subregionIndex,
                                       int basisChunk) const
{
  verifySubregionIndex(subregionIndex);
  TEST_FOR_EXCEPTION(basisChunk < 0 
                     || basisChunk >= nBasisChunks(subregionIndex), 
                     RuntimeError, 
                     "Basis chunk index=" << basisChunk
                     << " out of range [0, " << nBasisChunks(subregionIndex)
                     << ") on subregion=" << subregionIndex);
}



void DOFMapBase::getDOFsForCell(int cellDim, int cellLID,
                                int funcID,
                                Array<int>& dofs) const
{
  TimeMonitor timer(dofLookupTimer());
  
  Array<int> tmp = tuple(cellLID);
  RefCountPtr<Array<int> > cellLIDs = rcp(&tmp, false);
  Array<RefCountPtr<Array<int> > > batches;
  Array<int> subregionIndices;

  getHomogeneousCellBatches(cellLIDs, batches, subregionIndices);
  
  Array<Array<int> > allDofs;
  Array<int> nNodes;
  getDOFsForCellBatch(cellDim, *batches[0], allDofs, nNodes);

  int chunkNumber = chunkForFuncID(subregionIndices[0], funcID);
  int funcIndex = indexForFuncID(subregionIndices[0], chunkNumber, funcID);
  dofs.resize(nNodes[chunkNumber]);
  for (int i=0; i<nNodes[chunkNumber]; i++)
    {
      dofs[i] = allDofs[chunkNumber][nNodes[chunkNumber]*funcIndex + i];
    }
}

Time& DOFMapBase::dofLookupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("unbatched dof lookup"); 
  return *rtn;
}

Time& DOFMapBase::batchedDofLookupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("batched dof lookup"); 
  return *rtn;
}

