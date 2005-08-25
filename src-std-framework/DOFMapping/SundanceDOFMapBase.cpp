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

#include "SundanceDOFMapBase.hpp"
#include "Teuchos_MPIContainerComm.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



DOFMapBase::DOFMapBase(const Mesh& mesh)
  : localProcID_(mesh.comm().getRank()),
    mesh_(mesh),
    cellSets_(),
    funcIDOnCellSets_(),
    cellDimOnCellSets_(),
    lowestLocalDOF_(),
    numDOFs_(),
    ghostIndices_(rcp(new Array<int>())),
    dofsHaveBeenAssigned_()
{;}

void DOFMapBase::getDOFsForCell(int cellDim, int cellLID,
                                int funcID,
                                Array<int>& dofs) const
{
  Array<int> allDofs;
  unsigned int nNodes;
  getDOFsForCellBatch(cellDim, tuple(cellLID), allDofs, nNodes);

  dofs.resize(nNodes);
  unsigned int nFuncs = allDofs.size()/nNodes;
  for (unsigned int i=0; i<nNodes; i++)
    {
      dofs[i] = allDofs[funcID + nFuncs*i];
    }
}
