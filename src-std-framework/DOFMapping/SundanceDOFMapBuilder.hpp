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

#ifndef SUNDANCE_DOFMAPBUILDER_H
#define SUNDANCE_DOFMAPBUILDER_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCFMeshPair.hpp"
#include "SundanceMap.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;

namespace Internal
{
using namespace Teuchos;

/** 
 * 
 */
class DOFMapBuilder : public TSFExtended::ParameterControlledObjectWithVerbosity<DOFMapBase>
{
public:
  /** */
  DOFMapBuilder(const ParameterList& verbParams=*DOFMapBase::defaultVerbParams());
  /** */
  DOFMapBuilder(const Mesh& mesh, const RefCountPtr<EquationSet>& eqn, 
    bool findBCCols, const ParameterList& verbParams);

  /** */
  const Array<RefCountPtr<DOFMapBase> >& rowMap() const {return rowMap_;}

  /** */
  const Array<RefCountPtr<DOFMapBase> >& colMap() const {return colMap_;}

  /** */
  const Array<RefCountPtr<Array<int> > >& isBCRow() const {return isBCRow_;}

  /** */
  const Array<RefCountPtr<Array<int> > >& isBCCol() const {return isBCCol_;}


  /** */
  const Array<RefCountPtr<std::set<int> > >& remoteBCCols() const 
    {return remoteBCCols_;}

  Array<Array<BasisFamily> > testBasisArray() const ;

  Array<Array<BasisFamily> > unkBasisArray() const ;

  Array<Array<Set<CellFilter> > > testCellFilters() const ;

  Array<Array<Set<CellFilter> > > unkCellFilters() const ;

  const Mesh& mesh() const {return mesh_;}



  RefCountPtr<DOFMapBase> makeMap(const Mesh& mesh,
    const Array<BasisFamily>& basis,
    const Array<Set<CellFilter> >& filters) ;

  bool hasOmnipresentNodalMap(const Array<BasisFamily>& basis,
    const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool hasCommonDomain(const Array<Set<CellFilter> >& filters) const ;

  bool hasHomogeneousBasis(const Array<BasisFamily>& basis) const ;

  bool hasNodalBasis(const Array<BasisFamily>& basis) const ;

  bool hasCellBasis(const Array<BasisFamily>& basis) const ;

  bool allFuncsAreOmnipresent(const Mesh& mesh,
    const Array<Set<CellFilter> >& filters) const ;

  bool isWholeDomain(const Mesh& mesh,
    const Set<CellFilter>& filters) const ;

  CellFilter getMaxCellFilter(const Array<Set<CellFilter> >& filters) const ;

  static bool& allowNodalMap() {static bool rtn=true; return rtn;}

  /** */
  void extractUnkSetsFromEqnSet(const EquationSet& eqn,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  void extractVarSetsFromEqnSet(const EquationSet& eqn,
    Array<Set<int> >& funcSets,
    Array<CellFilter>& regions) const ;

  /** */
  SundanceUtils::Map<Set<int>, Set<CellFilter> > 
  buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
    const Array<CellFilter>& regions,
    const Mesh& mesh) const ;
        
  void getSubdomainUnkFuncMatches(const EquationSet& eqn,
    Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap) const ;
        
  void getSubdomainVarFuncMatches(const EquationSet& eqn,
    Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap) const ;

  Array<SundanceUtils::Map<Set<int>, CellFilter> > 
  funcDomains(const Mesh& mesh,
    const SundanceUtils::Map<CellFilter, Set<int> >& fmap,
    SundanceUtils::Map<CellFilter, SundanceUtils::Map<Set<int>, CellSet> >& inputToChildrenMap) const ;

  SundanceUtils::Map<CellFilter, Set<int> > domainToFuncSetMap(const Array<Set<CellFilter> >& filters) const ;

private:

  Set<CellFilter> reduceCellFilters(const Mesh& mesh,
    const Set<CellFilter>& inputSet) const ;

  bool hasUnks() const ;

  bool unksAreOmnipresent() const ;

  bool testsAreOmnipresent() const ;

  bool regionIsMaximal(int r) const ;

  bool isSymmetric(int block) const ;

  void markBCRows(int block) ;

  void markBCCols(int block) ;

  const MPIComm& comm() const {return mesh().comm();}

  void init(bool findBCCols);

  Mesh mesh_;

  RefCountPtr<EquationSet> eqn_;

  Array<RefCountPtr<DOFMapBase> > rowMap_;

  Array<RefCountPtr<DOFMapBase> > colMap_;

  Array<RefCountPtr<Array<int> > > isBCRow_;

  Array<RefCountPtr<Array<int> > > isBCCol_;

  Array<RefCountPtr<std::set<int> > > remoteBCCols_;
};
}
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
