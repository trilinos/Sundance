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
#include "SundanceLagrange.hpp"
#include "SundanceMixedDOFMap.hpp"
#include "SundanceNodalDOFMap.hpp"
#include "SundancePartialElementDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceInhomogeneousNodalDOFMap.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCFMeshPair.hpp"
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

static Time& cellFilterReductionTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("cell filter reduction"); 
  return *rtn;
}

static Time& findFuncDomainTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("finding func domains"); 
  return *rtn;
}

DOFMapBuilder::DOFMapBuilder(const Mesh& mesh, 
  const RefCountPtr<EquationSet>& eqn, bool findBCCols,
  const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<DOFMapBase>("DOF Map", verbParams),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    isBCRow_(),
    isBCCol_(),
    remoteBCCols_()
{
  init(findBCCols);
}

DOFMapBuilder::DOFMapBuilder()
  : mesh_(),
    eqn_(),
    rowMap_(),
    colMap_(),
    isBCRow_(),
    isBCCol_(),
    remoteBCCols_()
{}

RefCountPtr<DOFMapBase> DOFMapBuilder::makeMap(const Mesh& mesh,
                                               const Array<BasisFamily>& basis,
                                               const Array<Set<CellFilter> >& filters) 
{
  TimeMonitor timer(DOFBuilderCtorTimer());
  RefCountPtr<DOFMapBase> rtn;

  if (allowNodalMap() && hasOmnipresentNodalMap(basis, mesh, filters))
    {
      CellFilter maxCells = getMaxCellFilter(filters);
      rtn = rcp(new NodalDOFMap(mesh, basis.size(), maxCells));
    }
  else if (hasCellBasis(basis) && hasCommonDomain(filters))
    {
      TEST_FOR_EXCEPTION(filters[0].size() != 1, RuntimeError,
                         "only a single domain expected in construction of an element "
                         "DOF map");
      rtn = rcp(new PartialElementDOFMap(mesh, *filters[0].begin(), basis.size()));
    }
  else if (allFuncsAreOmnipresent(mesh, filters))
    {
      CellFilter maxCells = getMaxCellFilter(filters);
      rtn = rcp(new MixedDOFMap(mesh, basis, maxCells));
    }
  else if (hasNodalBasis(basis))
    {
      SundanceUtils::Map<CellFilter, Set<int> > fmap = domainToFuncSetMap(filters);
      SundanceUtils::Map<CellFilter, SundanceUtils::Map<Set<int>, CellSet> > inputChildren;

      Array<SundanceUtils::Map<Set<int>, CellFilter> > disjoint 
        = DOFMapBuilder::funcDomains(mesh, fmap, inputChildren);

      rtn = rcp(new InhomogeneousNodalDOFMap(mesh, disjoint));
    }
  else
    {
      TEST_FOR_EXCEPT(true);
    }
  return rtn;
}


SundanceUtils::Map<CellFilter, Set<int> > DOFMapBuilder::domainToFuncSetMap(const Array<Set<CellFilter> >& filters) 
{
  Map<CellFilter, Set<int> > rtn;
  for (unsigned int i=0; i<filters.size(); i++)
    {
      const Set<CellFilter>& s = filters[i];
      for (Set<CellFilter>::const_iterator j=s.begin(); j!=s.end(); j++)
        {
          const CellFilter& cf = *j;
          if (rtn.containsKey(cf)) 
            {
              rtn[cf].put(i);
            }
          else
            {
              
              rtn.put(cf, makeSet((int) i));
            }
        }
    }
  return rtn;
}


void DOFMapBuilder
::getSubdomainUnkFuncMatches(const EquationSet& eqn,
                             Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap)
{
  fmap.resize(eqn.numUnkBlocks());
  
  for (unsigned int r=0; r<eqn.numRegions(); r++)
    {
      CellFilter subreg = eqn.region(r);
      Set<int> funcs = eqn.unksOnRegion(r).setUnion(eqn.bcUnksOnRegion(r));
      for (Set<int>::const_iterator i=funcs.begin(); i!=funcs.end(); i++)
        {
          int block = eqn.blockForUnkID(*i);
          if (fmap[block].containsKey(subreg))
            {
              fmap[block][subreg].put(*i);
            }
          else
            {
              fmap[block].put(subreg, makeSet(*i));
            }
        }
    }
}

void DOFMapBuilder
::getSubdomainVarFuncMatches(const EquationSet& eqn,
                             Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap)
{
  fmap.resize(eqn.numVarBlocks());
  
  for (unsigned int r=0; r<eqn.numRegions(); r++)
    {
      CellFilter subreg = eqn.region(r);
      Set<int> funcs = eqn.varsOnRegion(r).setUnion(eqn.bcVarsOnRegion(r));
      for (Set<int>::const_iterator i=funcs.begin(); i!=funcs.end(); i++)
        {
          int block = eqn.blockForVarID(*i);
          if (fmap[block].containsKey(subreg))
            {
              fmap[block][subreg].put(*i);
            }
          else
            {
              fmap[block].put(subreg, makeSet(*i));
            }
        }
    }
}

Array<SundanceUtils::Map<Set<int>, CellFilter> > DOFMapBuilder
::funcDomains(const Mesh& mesh,
              const SundanceUtils::Map<CellFilter, Set<int> >& fmap,
              SundanceUtils::Map<CellFilter, SundanceUtils::Map<Set<int>, CellSet> >& inputToChildrenMap)
{
  TimeMonitor timer(findFuncDomainTimer());
  Array<Array<CellFilter> > filters(mesh.spatialDim()+1);
  Array<Array<Set<int> > > funcs(mesh.spatialDim()+1);

  for (SundanceUtils::Map<CellFilter, Set<int> >::const_iterator 
         i=fmap.begin(); i!=fmap.end(); i++)
    {
      int d = i->first.dimension(mesh);
      filters[d].append(i->first);
      funcs[d].append(i->second);
    }
  Array<Array<CFMeshPair> > tmp(mesh.spatialDim()+1);
  for (unsigned int d=0; d<tmp.size(); d++)
    {
      if (filters[d].size() != 0U)
        tmp[d] = findDisjointFilters(filters[d], funcs[d], mesh);
    }

  for (unsigned int d=0; d<tmp.size(); d++)
    {
      for (unsigned int r=0; r<tmp[d].size(); r++)
        {
          for (unsigned int p=0; p<filters[d].size(); p++)
            {
              if (tmp[d][r].filter().isSubsetOf(filters[d][p], mesh)) 
                {
                  if (inputToChildrenMap.containsKey(filters[d][p]))
                    {
                      SundanceUtils::Map<Set<int>, CellSet>& m 
                        = inputToChildrenMap[filters[d][p]];
                      if (m.containsKey(tmp[d][r].funcs()))
                        {
                          m.put(tmp[d][r].funcs(), m[tmp[d][r].funcs()].setUnion(tmp[d][r].cellSet())); 
                        }
                      else
                        {
                          m.put(tmp[d][r].funcs(), tmp[d][r].cellSet()); 
                        }
                    }
                  else
                    {
                      SundanceUtils::Map<Set<int>, CellSet> m;
                      m.put(tmp[d][r].funcs(), tmp[d][r].cellSet());
                      inputToChildrenMap.put(filters[d][p], m);
                    }
                }
            }
        }
    }

  Array<SundanceUtils::Map<Set<int>, CellFilter> > rtn(mesh.spatialDim()+1);
  for (unsigned int d=0; d<tmp.size(); d++)
    {
      if (tmp[d].size() == 0U) continue;
      for (unsigned int i=0; i<tmp[d].size(); i++)
        {
          const Set<int>& f = tmp[d][i].funcs();
          const CellFilter& cf = tmp[d][i].filter();
          if (rtn[d].containsKey(f))
            {
              rtn[d].put(f, rtn[d][f] + cf);
            }
          else
            {
              rtn[d].put(f, cf);
            }
        }
    }

  return rtn;
}


void DOFMapBuilder::init(bool findBCCols)
{
  rowMap_.resize(eqn_->numVarBlocks());
  colMap_.resize(eqn_->numUnkBlocks());
  isBCRow_.resize(eqn_->numVarBlocks());
  isBCCol_.resize(eqn_->numUnkBlocks());

  Array<Array<BasisFamily> > testBasis = testBasisArray();
  Array<Array<Set<CellFilter> > > testRegions = testCellFilters();

  for (unsigned int br=0; br<eqn_->numVarBlocks(); br++)
    {
      rowMap_[br] = makeMap(mesh_, testBasis[br], testRegions[br]);
      markBCRows(br);
    }      


  Array<Array<BasisFamily> > unkBasis = unkBasisArray();
  Array<Array<Set<CellFilter> > > unkRegions = unkCellFilters();

  for (unsigned int bc=0; bc<eqn_->numUnkBlocks(); bc++)
    {
      if (isSymmetric(bc))
        {
          colMap_[bc] = rowMap_[bc];
        }
      else
        {
          colMap_[bc] = makeMap(mesh_, unkBasis[bc], unkRegions[bc]);
        }
      if (findBCCols) markBCCols(bc);
    }
}

void DOFMapBuilder::extractUnkSetsFromEqnSet(const EquationSet& eqn,
                                             Array<Set<int> >& funcSets,
                                             Array<CellFilter>& regions)
{
  funcSets.resize(eqn.numRegions());
  regions.resize(eqn.numRegions());
  for (unsigned int r=0; r<eqn.numRegions(); r++)
    {
      regions[r] = eqn.region(r);
      funcSets[r] = eqn.unksOnRegion(r).setUnion(eqn.bcUnksOnRegion(r));
    }
}

void DOFMapBuilder::extractVarSetsFromEqnSet(const EquationSet& eqn,
                                             Array<Set<int> >& funcSets,
                                             Array<CellFilter>& regions)
{
  funcSets.resize(eqn.numRegions());
  regions.resize(eqn.numRegions());
  for (unsigned int r=0; r<eqn.numRegions(); r++)
    {
      regions[r] = eqn.region(r);
      funcSets[r] = eqn.varsOnRegion(r).setUnion(eqn.bcVarsOnRegion(r));
    }
}

SundanceUtils::Map<Set<int>, Set<CellFilter> > 
DOFMapBuilder::buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
                                      const Array<CellFilter>& regions,
                                      const Mesh& mesh)
{
  SundanceUtils::Map<Set<int>, Set<CellFilter> > tmp;
  
  for (unsigned int r=0; r<regions.size(); r++)
    {
      const CellFilter& reg = regions[r];
      if (!tmp.containsKey(funcSets[r]))
        {
          tmp.put(funcSets[r], Set<CellFilter>());
        }
      tmp[funcSets[r]].put(reg);
    }
  
  /* eliminate overlap between cell filters */
  SundanceUtils::Map<Set<int>, Set<CellFilter> > rtn;
  for (SundanceUtils::Map<Set<int>, Set<CellFilter> >::const_iterator 
         i=tmp.begin(); i!=tmp.end(); i++)
    {
      rtn.put(i->first, reduceCellFilters(mesh, i->second));
    }

  return rtn;
}

bool DOFMapBuilder::hasOmnipresentNodalMap(const Array<BasisFamily>& basis,
                                           const Mesh& mesh,
                                           const Array<Set<CellFilter> >& filters)
{
  return hasNodalBasis(basis) && allFuncsAreOmnipresent(mesh, filters);
}
                                           
bool DOFMapBuilder::hasHomogeneousBasis(const Array<BasisFamily>& basis) 
{
  for (unsigned int i=1; i<basis.size(); i++)
    {
      if (!(basis[i] == basis[i-1])) return false;
    }
  return true;
}

                                           
bool DOFMapBuilder::hasCommonDomain(const Array<Set<CellFilter> >& filters) 
{
  Set<CellFilter> first = filters[0];
  for (unsigned int i=1; i<filters.size(); i++) 
    {
      if (! (filters[i] == first) ) return false;
    }
  return true;
}                           

bool DOFMapBuilder::hasNodalBasis(const Array<BasisFamily>& basis)
{
  for (unsigned int i=0; i<basis.size(); i++)
    {
      const Lagrange* lagr 
        = dynamic_cast<const Lagrange*>(basis[0].ptr().get());
      if (lagr==0 || basis[0].order()!=1) return false;
    }
  return true;
}


bool DOFMapBuilder::hasCellBasis(const Array<BasisFamily>& basis)
{
  for (unsigned int i=0; i<basis.size(); i++)
    {
      const Lagrange* lagr 
        = dynamic_cast<const Lagrange*>(basis[0].ptr().get());
      if (lagr==0 || basis[0].order()!=0) return false;
    }
  return true;
}

bool DOFMapBuilder::allFuncsAreOmnipresent(const Mesh& mesh, 
                                           const Array<Set<CellFilter> >& filters) 
{
  Set<Set<CellFilter> > distinctSets;
  for (unsigned int i=0; i<filters.size(); i++)
    {
      distinctSets.put(filters[i]);

    }
  for (Set<Set<CellFilter> >::const_iterator 
         iter=distinctSets.begin(); iter != distinctSets.end(); iter++)
    {
      if (!isWholeDomain(mesh, *iter)) return false;
    }

  return true;
}

bool DOFMapBuilder::isWholeDomain(const Mesh& mesh, 
                                  const Set<CellFilter>& filters) 
{
  CellFilter allMax = new MaximalCellFilter();
  CellSet remainder = allMax.getCells(mesh);

  for (Set<CellFilter>::const_iterator 
         i=filters.begin(); i!=filters.end(); i++)
    {
      const CellFilter& cf = *i;
      if (0 != dynamic_cast<const MaximalCellFilter*>(cf.ptr().get()))
        {
          return true;
        }
      if (cf.dimension(mesh) != mesh.spatialDim()) continue;
      CellSet cells = cf.getCells(mesh);
      remainder = remainder.setDifference(cells);
      if (remainder.begin() == remainder.end()) return true;
    }

  return false;
}



CellFilter DOFMapBuilder::getMaxCellFilter(const Array<Set<CellFilter> >& filters) 
{
  for (unsigned int i=0; i<filters.size(); i++)
    {
      const Set<CellFilter>& cfs = filters[i];
      if (cfs.size() != 1U) continue;
      const CellFilter& cf = *cfs.begin();
      if (0 != dynamic_cast<const MaximalCellFilter*>(cf.ptr().get()))
        return cf;
    }
  TEST_FOR_EXCEPT(true);
  return new MaximalCellFilter();
}



Array<Array<Set<CellFilter> > > DOFMapBuilder::testCellFilters() const
{
  Array<Array<Set<CellFilter> > > rtn(eqn_->numVarBlocks());

  for (unsigned int b=0; b<eqn_->numVarBlocks(); b++)
    {
      for (unsigned int i=0; i<eqn_->numVars(b); i++) 
        {
          int testID = eqn_->unreducedVarID(b, i);
          Set<CellFilter> s;
          const Set<OrderedHandle<CellFilterStub> >& cfs 
            = eqn_->regionsForTestFunc(testID);
          for (Set<OrderedHandle<CellFilterStub> >::const_iterator 
                 j=cfs.begin(); j!=cfs.end(); j++)
            {
              RefCountPtr<CellFilterBase> cfb 
                = rcp_dynamic_cast<CellFilterBase>(j->ptr());
              TEST_FOR_EXCEPT(cfb.get()==0);
              CellFilter cf = j->ptr();
              s.put(cf);
            }
          Set<CellFilter> reducedS = reduceCellFilters(mesh(), s);
          rtn[b].append(reducedS);
        }
    }
  return rtn;
}

Array<Array<Set<CellFilter> > > DOFMapBuilder::unkCellFilters() const
{
  Array<Array<Set<CellFilter> > > rtn(eqn_->numUnkBlocks());

  for (unsigned int b=0; b<eqn_->numUnkBlocks(); b++)
    {
      for (unsigned int i=0; i<eqn_->numUnks(b); i++) 
        {
          int unkID = eqn_->unreducedUnkID(b, i);
          Set<CellFilter> s;
          const Set<OrderedHandle<CellFilterStub> >& cfs 
            = eqn_->regionsForUnkFunc(unkID);
          for (Set<OrderedHandle<CellFilterStub> >::const_iterator 
                 j=cfs.begin(); j!=cfs.end(); j++)
            {
              RefCountPtr<CellFilterBase> cfb 
                = rcp_dynamic_cast<CellFilterBase>(j->ptr());
              TEST_FOR_EXCEPT(cfb.get()==0);
              CellFilter cf = j->ptr();
              s.put(cf);
            }
          Set<CellFilter> reducedS = reduceCellFilters(mesh(), s);
          rtn[b].append(reducedS);
        }
    }
  return rtn;
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


Set<CellFilter> DOFMapBuilder
::reduceCellFilters(const Mesh& mesh, 
                    const Set<CellFilter>& inputSet)  
{
  TimeMonitor timer(cellFilterReductionTimer());
  Set<CellFilter> rtn;
  /* If the input set explicitly contains all maximal cells, we're done */
  CellFilter m = new MaximalCellFilter();
  if (inputSet.contains(m))
    {
      rtn.put(m);
      return rtn;
    }

  /* Next, see if combining the maximal-dimension filters in the
   * input set gives us all maximal cells. */
  CellFilter myMaxFilters;
  for (Set<CellFilter>::const_iterator 
         i=inputSet.begin(); i!=inputSet.end(); i++)
    {
      CellFilter f = *i;
      if (f.dimension(mesh) != mesh.spatialDim()) continue;
      myMaxFilters = myMaxFilters + f;
    }
  CellSet allMax = m.getCells(mesh);
  CellSet myMax = myMaxFilters.getCells(mesh);

  CellSet diff = allMax.setDifference(myMax);
  /* if the difference between the collected max cell set and the known
   * set of all max cells is empty, then we're done */
  if (diff.begin() == diff.end())
    {
      rtn.put(m);
      return rtn;
    }
  
  /* Otherwise, we return the remaining max cell filters, and possibly
   * some lower-dimensional filters to be identified below. */
  if (myMax.begin() != myMax.end()) rtn.put(myMaxFilters);

  /* At this point, we can eliminate as redundant any lower-dimensional
   * cell filters all of whose cells are facets of our subset of
   * maximal cell filters. Any other lower-dim cell filters must be 
   * appended to our list. */
  for (Set<CellFilter>::const_iterator 
         i=inputSet.begin(); i!=inputSet.end(); i++)
    {
      CellFilter f = *i;
      if (f.dimension(mesh) == mesh.spatialDim()) continue;
      CellSet s = f.getCells(mesh);
      if (s.areFacetsOf(myMax)) continue;
      /* if we're here, then we have a lower-dimensional cell filter
       * whose cells are not facets of cells in our maximal cell filters.
       * These will need to be processed separately in assigning DOFs.
       */
      rtn.put(f);
    }
  return rtn;
}


bool DOFMapBuilder::isSymmetric(int b) const 
{
  if ((int)eqn_->numVarBlocks() < b || (int)eqn_->numUnkBlocks() < b) return false;

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

  RefCountPtr<Array<int> > cellLID = rcp(new Array<int>());
  Array<RefCountPtr<Array<int> > > cellBatches;
  const RefCountPtr<DOFMapBase>& rowMap = rowMap_[block];

  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      /* find the cells in this region */
      CellFilter region = eqn_->region(r);

      if (!eqn_->isBCRegion(r)) continue;

      int dim = region.dimension(mesh_);
      CellSet cells = region.getCells(mesh_);
      cellLID->resize(0);
      for (CellIterator c=cells.begin(); c != cells.end(); c++)
        {
          cellLID->append(*c);
        }
      if (cellLID->size() == 0U) continue;
      
      /* find the functions that appear in BCs on this region */
      const Set<int>& allBcFuncs = eqn_->bcVarsOnRegion(r);

      Set<int> bcFuncs;
      for (Set<int>::const_iterator 
             i=allBcFuncs.begin(); i != allBcFuncs.end(); i++)
        {
          if (block == eqn_->blockForVarID(*i)) 
            {
              bcFuncs.put(eqn_->reducedVarID(*i));
            }
        }
      if (bcFuncs.size()==0) continue;
      Array<int> bcFuncID = bcFuncs.elements();

      Array<Array<int> > dofs;
      Array<int> nNodes;

      RefCountPtr<const MapStructure> s 
        = rowMap->getDOFsForCellBatch(dim, *cellLID, bcFuncs, dofs, nNodes);
      int offset = rowMap->lowestLocalDOF();
      int high = offset + rowMap->numLocalDOFs();
      
      for (unsigned int c=0; c<cellLID->size(); c++)
        {
          for (int b=0; b< s->numBasisChunks(); b++)
            {
              int nFuncs = s->numFuncs(b);
              for (int n=0; n<nNodes[b]; n++)
                {
                  for (unsigned int f=0; f<bcFuncID.size(); f++)
                    {
                      int chunk = s->chunkForFuncID(bcFuncID[f]);
                      if (chunk != b) continue;
                      int funcOffset = s->indexForFuncID(bcFuncID[f]);
                      int dof = dofs[b][(c*nFuncs + funcOffset)*nNodes[b]+n];
                      if (dof < offset || dof >= high) continue;
                      (*isBCRow_[block])[dof-offset]=true;
                    }
                }
            }
        }
    }
}
        


void DOFMapBuilder::markBCCols(int block)
{
  isBCCol_[block] = rcp(new Array<int>(colMap_[block]->numLocalDOFs()));
  int ndof = colMap_[block]->numLocalDOFs();
  Array<int>& isBC = *isBCCol_[block];
  for (int i=0; i<ndof; i++) isBC[i] = false;

  RefCountPtr<Array<int> > cellLID = rcp(new Array<int>());
  Array<RefCountPtr<Array<int> > > cellBatches;
  const RefCountPtr<DOFMapBase>& colMap = colMap_[block];

  for (unsigned int r=0; r<eqn_->numRegions(); r++)
    {
      /* find the cells in this region */
      CellFilter region = eqn_->region(r);

      if (!eqn_->isBCRegion(r)) continue;

      int dim = region.dimension(mesh_);
      CellSet cells = region.getCells(mesh_);
      cellLID->resize(0);
      for (CellIterator c=cells.begin(); c != cells.end(); c++)
        {
          cellLID->append(*c);
        }
      if (cellLID->size() == 0U) continue;
      
      /* find the functions that appear in BCs on this region */
      const Set<int>& allBcFuncs = eqn_->bcUnksOnRegion(r);

      Set<int> bcFuncs;
      for (Set<int>::const_iterator 
             i=allBcFuncs.begin(); i != allBcFuncs.end(); i++)
        {
          if (block == eqn_->blockForUnkID(*i)) 
            {
              bcFuncs.put(eqn_->reducedUnkID(*i));
            }
        }
      if (bcFuncs.size()==0) continue;
      Array<int> bcFuncID = bcFuncs.elements();

      Array<Array<int> > dofs;
      Array<int> nNodes;

      RefCountPtr<const MapStructure> s 
        = colMap->getDOFsForCellBatch(dim, *cellLID, bcFuncs, dofs, nNodes);
      int offset = colMap->lowestLocalDOF();
      int high = offset + colMap->numLocalDOFs();
      
      for (unsigned int c=0; c<cellLID->size(); c++)
        {
          for (int b=0; b< s->numBasisChunks(); b++)
            {
              int nFuncs = s->numFuncs(b);
              for (int n=0; n<nNodes[b]; n++)
                {
                  for (unsigned int f=0; f<bcFuncID.size(); f++)
                    {
                      int chunk = s->chunkForFuncID(bcFuncID[f]);
                      if (chunk != b) continue;
                      int funcOffset = s->indexForFuncID(bcFuncID[f]);
                      int dof = dofs[b][(c*nFuncs + funcOffset)*nNodes[b]+n];
                      if (dof < offset || dof >= high) 
                      {
                        remoteBCCols_[block]->insert(dof);
                      }
                      else
                      {
                        (*isBCCol_[block])[dof-offset]=true;
                      }
                    }
                }
            }
        }
    }
}
        
