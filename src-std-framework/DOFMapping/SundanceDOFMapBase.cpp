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
#include "SundanceDOFMapBase.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;


static Time& dofLookupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("unbatched dof lookup"); 
  return *rtn;
}

DOFMapBase::DOFMapBase(const Mesh& mesh,
                       const Array<BasisFamily>& basis)
  : localProcID_(mesh.comm().getRank()),
    mesh_(mesh),
    cellSets_(),
    funcIDOnCellSets_(),
    cellDimOnCellSets_(),
    lowestLocalDOF_(),
    numDOFs_(),
    ghostIndices_(rcp(new Array<int>())),
    dofsHaveBeenAssigned_(),
    chunkBasis_(),
    chunkFuncIDs_(),
    funcIDToChunkMap_(basis.size()),
    funcIDToIndexMap_(basis.size())
{
  SundanceUtils::Map<BasisFamily, int> basisToChunkMap;
  
  int nBasis = basis.size();
  int chunk = 0;
  for (int i=0; i<nBasis; i++)
    {
      if (!basisToChunkMap.containsKey(basis[i]))
        {
          chunkBasis_.append(basis[i]);
          basisToChunkMap.put(basis[i], chunk);
          chunkFuncIDs_.append(tuple(i));
          chunk++;
        }
      else
        {
          int b = basisToChunkMap.get(basis[i]);
          chunkFuncIDs_[b].append(i);
        }

      funcIDToChunkMap_[i] = basisToChunkMap.get(basis[i]);
      funcIDToIndexMap_[i] = chunkFuncIDs_[funcIDToChunkMap_[i]].size()-1;
    }

  /* identify all functions as existing on the maximal cell set */
  Array<int> fid(basis.size());
  for (int f=0; f<nBasis; f++) fid[f] = f;
  funcIDOnCellSets().append(fid);
}

void DOFMapBase::getDOFsForCell(int cellDim, int cellLID,
                                int funcID,
                                Array<int>& dofs) const
{
  TimeMonitor timer(dofLookupTimer());
  
  Array<Array<int> > allDofs;
  Array<int> nNodes;
  getDOFsForCellBatch(cellDim, tuple(cellLID), allDofs, nNodes);

  int chunkNumber = chunkForFuncID(funcID);
  int funcIndex = indexForFuncID(funcID);
  dofs.resize(nNodes[chunkNumber]);
  for (int i=0; i<nNodes[chunkNumber]; i++)
    {
      dofs[i] = allDofs[chunkNumber][nNodes[chunkNumber]*funcIndex + i];
    }
}
