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
#include "SundanceHomogeneousDOFMap.hpp"
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
    isBCRow_(rcp(new Array<int>()))
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
  /* If the lists of test and unknown functions have the same basis families
   * in the same order, then the row and column spaces can share a single DOF
   * map. */
  if (isSymmetric())
    {
      /* if every test function is defined on the maximal cell set, and if
       * all test functions share a common basis, then we 
       * can build a homogeneous DOF map. */
      if (testsAreHomogeneous() && testsAreOmnipresent())
        {
          Expr test0 = eqn_->varFunc(0);
          BasisFamily basis0 = BasisFamily::getBasis(test0);
          rowMap_ = rcp(new HomogeneousDOFMap(mesh_, basis0, 
                                              eqn_->numVars()));
          colMap_ = rowMap_;
        }
      else
        {
          SUNDANCE_ERROR("DOFMapBuilder::init() non-homogeneous test function spaces not yet supported");
        }
    }
  else
    {
      /* if every test function is defined on the maximal cell set, and if
       * all test functions share a common basis, then we 
       * can build a homogeneous DOF map. */
      if (testsAreHomogeneous() && testsAreOmnipresent())
        {
          Expr test0 = eqn_->varFunc(0);
          BasisFamily basis0 = BasisFamily::getBasis(test0);
          rowMap_ = rcp(new HomogeneousDOFMap(mesh_, basis0,  
                                              eqn_->numVars()));
        }
      else
        {
          SUNDANCE_ERROR("DOFMapBuilder::init() non-homogeneous test function spaces not yet supported");
        }
      
      if (hasUnks())
        {
          /* if every unk function is defined on the maximal cell set, and if
           * all unk functions share a common basis, then we 
           * can build a homogeneous DOF map. */
          if (unksAreHomogeneous() && unksAreOmnipresent())
            {
              Expr unk0 = eqn_->unkFunc(0);
              BasisFamily basis0 = BasisFamily::getBasis(unk0);
              colMap_ = rcp(new HomogeneousDOFMap(mesh_, basis0, 
                                                  eqn_->numUnks()));
            }
          else
            {
              SUNDANCE_ERROR("DOFMapBuilder::init() non-homogeneous unknown function spaces not yet supported");
            }
        }
    }
  markBCRows();
}

Array<BasisFamily> DOFMapBuilder::testBasisArray() const 
{
  Array<BasisFamily> rtn;
  for (unsigned int i=0; i<eqn_->numVars(); i++) 
    {
      rtn.append(BasisFamily::getBasis(eqn_->varFunc(i)));
    }
  return rtn;
}

Array<BasisFamily> DOFMapBuilder::unkBasisArray() const 
{
  Array<BasisFamily> rtn;
  for (unsigned int i=0; i<eqn_->numUnks(); i++) 
    {
      rtn.append(BasisFamily::getBasis(eqn_->unkFunc(i)));
    }
  return rtn;
}

bool DOFMapBuilder::hasUnks() const
{
  return eqn_->numUnks() > 0;
}

bool DOFMapBuilder::unksAreHomogeneous() const 
{
  if (eqn_->numUnks() > 1)
    {
      BasisFamily basis0 = BasisFamily::getBasis(eqn_->unkFunc(0));
      for (unsigned int i=1; i<eqn_->numUnks(); i++) 
        {
          BasisFamily basis = BasisFamily::getBasis(eqn_->unkFunc(i));
          if (!(basis == basis0)) return false;
        }
    }
  return true;
}

bool DOFMapBuilder::testsAreHomogeneous() const 
{
  if (eqn_->numVars() > 1)
    {
      BasisFamily basis0 = BasisFamily::getBasis(eqn_->varFunc(0));
      for (unsigned int i=1; i<eqn_->numVars(); i++) 
        {
          BasisFamily basis = BasisFamily::getBasis(eqn_->varFunc(i));
          if (!(basis == basis0)) return false;
        }
    }
  return true;
}

bool DOFMapBuilder::unksAreOmnipresent() const
{
  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      if (regionIsMaximal(r))
        {
          if (eqn_->unksOnRegion(r).size() == eqn_->numUnks()) return true;
          else return false;
        }
    }
  return false;
}

bool DOFMapBuilder::testsAreOmnipresent() const
{
  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      if (regionIsMaximal(r))
        {
          if (eqn_->varsOnRegion(r).size() == eqn_->numVars()) return true;
          else return false;
        }
    }
  return false;
}


bool DOFMapBuilder::isSymmetric() const 
{
  if (eqn_->numVars() != eqn_->numUnks()) return false;

  for (unsigned int i=0; i<eqn_->numVars(); i++) 
    {
      BasisFamily basis1 = BasisFamily::getBasis(eqn_->varFunc(i));
      BasisFamily basis2 = BasisFamily::getBasis(eqn_->unkFunc(i));
      if (!(basis1 == basis2)) return false;
    }
  return true;
}

bool DOFMapBuilder::regionIsMaximal(int r) const 
{
  const CellFilterStub* reg = eqn_->region(r).get();
  return (dynamic_cast<const MaximalCellFilter*>(reg) != 0);
}

void DOFMapBuilder::markBCRows()
{
  isBCRow_->resize(rowMap_->numLocalDOFs());
  int ndof = rowMap_->numLocalDOFs();
  Array<int>& isBC = *isBCRow_;
  for (int i=0; i<ndof; i++) isBC[i] = false;
  Array<int> dofs;
  Array<int> cellLID;

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
      unsigned int nTestNodes;
      /* find the functions that appear in BCs on this region */
      const Set<int>& bcFuncs = eqn_->bcVarsOnRegion(r);
      Array<int> bcFuncID = bcFuncs.elements();
      for (unsigned int f=0; f<bcFuncID.size(); f++) 
        {
          bcFuncID[f] = eqn_->reducedVarID(bcFuncID[f]);
        }

      rowMap_->getDOFsForCellBatch(dim, cellLID, bcFuncID, dofs, nTestNodes);
      int offset = rowMap_->lowestLocalDOF();
      int high = offset + rowMap_->numLocalDOFs();
      for (unsigned int n=0; n<dofs.size(); n++) 
        {
          if (dofs[n] < offset || dofs[n] >= high) continue;
          (*isBCRow_)[dofs[n]-offset]=true;
        }
    }
}
