/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;


DOFMapBuilder::DOFMapBuilder(const Mesh& mesh, 
                             const RefCountPtr<EquationSet>& eqn)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    bcRows_(rcp(new Set<int>()))
{
  init();
}

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
          Expr test0 = eqn_->testFunc(0);
          BasisFamily basis0 = BasisFamily::getBasis(test0);
          rowMap_ = rcp(new HomogeneousDOFMap(mesh_, basis0, 
                                              eqn_->numTests()));
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
          Expr test0 = eqn_->testFunc(0);
          BasisFamily basis0 = BasisFamily::getBasis(test0);
          rowMap_ = rcp(new HomogeneousDOFMap(mesh_, basis0,  
                                              eqn_->numTests()));
        }
      else
        {
          SUNDANCE_ERROR("DOFMapBuilder::init() non-homogeneous test function spaces not yet supported");
        }

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
  markBCRows();
}
  
bool DOFMapBuilder::unksAreHomogeneous() const 
{
  if (eqn_->numUnks() > 1)
    {
      BasisFamily basis0 = BasisFamily::getBasis(eqn_->unkFunc(0));
      for (int i=1; i<eqn_->numUnks(); i++) 
        {
          BasisFamily basis = BasisFamily::getBasis(eqn_->unkFunc(i));
          if (!(basis == basis0)) return false;
        }
    }
  return true;
}

bool DOFMapBuilder::testsAreHomogeneous() const 
{
  if (eqn_->numTests() > 1)
    {
      BasisFamily basis0 = BasisFamily::getBasis(eqn_->testFunc(0));
      for (int i=1; i<eqn_->numTests(); i++) 
        {
          BasisFamily basis = BasisFamily::getBasis(eqn_->testFunc(i));
          if (!(basis == basis0)) return false;
        }
    }
  return true;
}

bool DOFMapBuilder::unksAreOmnipresent() const
{
  for (int r=0; r<eqn_->numRegions(); r++)
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
  for (int r=0; r<eqn_->numRegions(); r++)
    {
      if (regionIsMaximal(r))
        {
          if (eqn_->testsOnRegion(r).size() == eqn_->numTests()) return true;
          else return false;
        }
    }
  return false;
}


bool DOFMapBuilder::isSymmetric() const 
{
  if (eqn_->numTests() != eqn_->numUnks()) return false;

  for (int i=0; i<eqn_->numTests(); i++) 
    {
      BasisFamily basis1 = BasisFamily::getBasis(eqn_->testFunc(i));
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
  for (int r=0; r<eqn_->numRegions(); r++)
    {
      if (!eqn_->isBCRegion(r)) continue;

      /* find the cells in this region */
      CellFilter region = eqn_->region(r);
      int dim = region.dimension(mesh_);
      CellSet cells = region.getCells(mesh_);
      
      cerr << "finding bc rows for cell filter " << region << endl;
      cerr << "dim is " << dim << endl;
      cerr << "cells are " << cells << endl;

      /* find the functions that appear in BCs on this region */
      const Set<int>& bcFuncs = eqn_->bcTestsOnRegion(r);

      for (CellIterator c = cells.begin(); c != cells.end(); c++)
        {
          int cellLID = *c;
          for (Set<int>::const_iterator f = bcFuncs.begin(); 
               f != bcFuncs.end(); f++)
            {
              Array<int> dofs;
              int fid = eqn_->reducedTestID(*f);
              rowMap_->getDOFsForCell(dim, cellLID, fid, dofs);
              cerr << "cell=" << cellLID << " dofs=" << dofs << endl;
              for (int i=0; i<dofs.size(); i++) bcRows_->put(dofs[i]);
            }
        }
    }
}
