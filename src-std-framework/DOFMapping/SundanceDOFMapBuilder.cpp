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

#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceMixedDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;

static Time& DOFBuilderCtorTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("DOF map building"); 
  return *rtn;
}

DOFMapBuilder::DOFMapBuilder(const Mesh& mesh, 
                             const RefCountPtr<EquationSet>& eqn)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    isBCRow_()
{
  TimeMonitor timer(DOFBuilderCtorTimer());
  init();
}

DOFMapBuilder::DOFMapBuilder()
  : mesh_(),
    eqn_(),
    rowMap_(),
    colMap_(),
    isBCRow_()
{}

void DOFMapBuilder::init()
{
  rowMap_.resize(eqn_->numVarBlocks());
  colMap_.resize(eqn_->numUnkBlocks());
  isBCRow_.resize(eqn_->numVarBlocks());

  Array<Array<BasisFamily> > testBasis = testBasisArray();
  for (unsigned int br=0; br<eqn_->numVarBlocks(); br++)
    {
      /* if every test function is defined on the maximal cell set, then we 
       * can build a mixed DOF map. */
      if (testsAreOmnipresent())
        {
          rowMap_[br] = rcp(new MixedDOFMap(mesh_, testBasis[br]));
        }
      else
        {
          SUNDANCE_ERROR("DOFMapBuilder::init() non-omnipresent test function "
                         "spaces not yet supported");
        }
      markBCRows(br);
    }


  Array<Array<BasisFamily> > unkBasis = unkBasisArray();
  for (unsigned int bc=0; bc<eqn_->numUnkBlocks(); bc++)
    {
      if (isSymmetric(bc))
        {
          colMap_[bc] = rowMap_[bc];
        }
      else
        {
          /* if every unk function is defined on the maximal cell set, then we 
           * can build a mixed DOF map. */
          if (unksAreOmnipresent())
            {
              colMap_[bc] = rcp(new MixedDOFMap(mesh_, unkBasis[bc]));
            }
          else
            {
              SUNDANCE_ERROR("DOFMapBuilder::init() non-omnipresent "
                             "test function "
                             "spaces not yet supported");
            }
        }
    }
}

Array<Array<BasisFamily> > DOFMapBuilder::testBasisArray() const 
{
  Array<Array<BasisFamily> > rtn(eqn_->numVarBlocks());
  for (unsigned int b=0; b<eqn_->numVarBlocks(); b++)
    {
      for (unsigned int i=0; i<eqn_->numVars(b); i++) 
        {
          rtn[b].append(BasisFamily::getBasis(eqn_->varFunc(b, i)));
        }
    }
  return rtn;
}

Array<Array<BasisFamily> > DOFMapBuilder::unkBasisArray() const 
{
  Array<Array<BasisFamily> > rtn(eqn_->numUnkBlocks());
  for (unsigned int b=0; b<eqn_->numUnkBlocks(); b++)
    {
      for (unsigned int i=0; i<eqn_->numUnks(b); i++) 
        {
          rtn[b].append(BasisFamily::getBasis(eqn_->unkFunc(b, i)));
        }
    }
  return rtn;
}

bool DOFMapBuilder::unksAreOmnipresent() const
{
  int numUnks = 0;
  for (unsigned int b=0; b<eqn_->numUnkBlocks(); b++) 
    {
      numUnks += eqn_->numUnks(b);
    }

  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      if (regionIsMaximal(r))
        {
          if (eqn_->unksOnRegion(r).size() == numUnks) return true;
          else return false;
        }
    }
  return false;
}

bool DOFMapBuilder::testsAreOmnipresent() const
{
  int numVars = 0;
  for (unsigned int b=0; b<eqn_->numVarBlocks(); b++) 
    {
      numVars += eqn_->numVars(b);
    }
  
  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      if (regionIsMaximal(r))
        {
          if (eqn_->varsOnRegion(r).size() == numVars) return true;
          else return false;
        }
    }
  return false;
}

bool DOFMapBuilder::isSymmetric(int b) const 
{
  if (eqn_->numVarBlocks() < b || eqn_->numUnkBlocks() < b) return false;

  if (eqn_->numVars(b) != eqn_->numUnks(b)) return false;

  for (unsigned int i=0; i<eqn_->numVars(b); i++) 
    {
      BasisFamily basis1 = BasisFamily::getBasis(eqn_->varFunc(b,i));
      BasisFamily basis2 = BasisFamily::getBasis(eqn_->unkFunc(b,i));
      if (!(basis1 == basis2)) return false;
    }
  return true;
}

bool DOFMapBuilder::regionIsMaximal(int r) const 
{
  const CellFilterStub* reg = eqn_->region(r).get();
  return (dynamic_cast<const MaximalCellFilter*>(reg) != 0);
}

void DOFMapBuilder::markBCRows(int block)
{
  isBCRow_[block] = rcp(new Array<int>(rowMap_[block]->numLocalDOFs()));
  int ndof = rowMap_[block]->numLocalDOFs();
  Array<int>& isBC = *isBCRow_[block];
  for (int i=0; i<ndof; i++) isBC[i] = false;
  Array<Array<int> > dofs(rowMap_[block]->nChunks());
  Array<int> cellLID;
  const RefCountPtr<DOFMapBase>& rowMap = rowMap_[block];

  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      /* find the cells in this region */
      CellFilter region = eqn_->region(r);

      if (!eqn_->isBCRegion(r)) continue;

      int dim = region.dimension(mesh_);
      CellSet cells = region.getCells(mesh_);
      cellLID.resize(0);
      for (CellIterator c=cells.begin(); c != cells.end(); c++)
        {
          cellLID.append(*c);
        }

      /* find the functions that appear in BCs on this region */
      const Set<int>& allBcFuncs = eqn_->bcVarsOnRegion(r);
      Set<int> bcFuncs;
      for (Set<int>::const_iterator 
             i=allBcFuncs.begin(); i != allBcFuncs.end(); i++)
        {
          if (block == eqn_->blockForVarID(*i)) bcFuncs.put(*i);
        }
      if (bcFuncs.size()==0) continue;
      Array<int> bcFuncID = bcFuncs.elements();
      for (unsigned int f=0; f<bcFuncID.size(); f++) 
        {
          bcFuncID[f] = eqn_->reducedVarID(bcFuncID[f]);
        }

      Array<int> nNodes;
      rowMap->getDOFsForCellBatch(dim, cellLID, dofs, nNodes);
      int offset = rowMap->lowestLocalDOF();
      int high = offset + rowMap->numLocalDOFs();
      if (cellLID.size()==0) continue;
      for (unsigned int c=0; c<cellLID.size(); c++)
        {
          for (int b=0; b<rowMap->nChunks(); b++)
            {
              int nFuncs = rowMap->nFuncs(b);
              for (int n=0; n<nNodes[b]; n++)
                {
                  for (unsigned int f=0; f<bcFuncID.size(); f++)
                    {
                      int chunk = rowMap->chunkForFuncID(bcFuncID[f]);
                      if (chunk != b) continue;
                      int funcOffset = rowMap->indexForFuncID(bcFuncID[f]);
                      int dof = dofs[b][(c*nFuncs + funcOffset)*nNodes[b]+n];
                      if (dof < offset || dof >= high) continue;
                      (*isBCRow_[block])[dof-offset]=true;
                    }
                }
            }
        }
    }
}
        
