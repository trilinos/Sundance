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
    class DOFMapBuilder : public TSFExtended::ObjectWithVerbosity<DOFMapBase>
    {
    public:
      /** */
      DOFMapBuilder();
      /** */
      DOFMapBuilder(const Mesh& mesh, const RefCountPtr<EquationSet>& eqn);

      /** */
      const Array<RefCountPtr<DOFMapBase> >& rowMap() const {return rowMap_;}

      /** */
      const Array<RefCountPtr<DOFMapBase> >& colMap() const {return colMap_;}

      /** */
      const Array<RefCountPtr<Array<int> > >& isBCRow() const {return isBCRow_;}

      Array<Array<BasisFamily> > testBasisArray() const ;

      Array<Array<BasisFamily> > unkBasisArray() const ;

      Array<Array<Set<CellFilter> > > testCellFilters() const ;

      Array<Array<Set<CellFilter> > > unkCellFilters() const ;

      const Mesh& mesh() const {return mesh_;}



      static RefCountPtr<DOFMapBase> makeMap(const Mesh& mesh,
                                               const Array<BasisFamily>& basis,
                                             const Array<Set<CellFilter> >& filters) ;

      static bool hasOmnipresentNodalMap(const Array<BasisFamily>& basis,
                                         const Mesh& mesh,
                                         const Array<Set<CellFilter> >& filters) ;

      static bool hasCommonDomain(const Array<Set<CellFilter> >& filters) ;

      static bool hasHomogeneousBasis(const Array<BasisFamily>& basis) ;

      static bool hasNodalBasis(const Array<BasisFamily>& basis) ;

      static bool hasCellBasis(const Array<BasisFamily>& basis) ;

      static bool allFuncsAreOmnipresent(const Mesh& mesh,
                                         const Array<Set<CellFilter> >& filters);

      static bool isWholeDomain(const Mesh& mesh,
                                const Set<CellFilter>& filters);

      static CellFilter getMaxCellFilter(const Array<Set<CellFilter> >& filters);

      static bool& allowNodalMap() {static bool rtn=true; return rtn;}

      /** */
      static void extractUnkSetsFromEqnSet(const EquationSet& eqn,
                                           Array<Set<int> >& funcSets,
                                           Array<CellFilter>& regions);

      /** */
      static void extractVarSetsFromEqnSet(const EquationSet& eqn,
                                           Array<Set<int> >& funcSets,
                                           Array<CellFilter>& regions);

      /** */
      static SundanceUtils::Map<Set<int>, Set<CellFilter> > 
      buildFuncSetToCFSetMap(const Array<Set<int> >& funcSets,
                             const Array<CellFilter>& regions,
                                      const Mesh& mesh);
        
      static void getSubdomainUnkFuncMatches(const EquationSet& eqn,
                                             Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap);
        
      static void getSubdomainVarFuncMatches(const EquationSet& eqn,
                                             Array<SundanceUtils::Map<CellFilter, Set<int> > >& fmap);

      static Array<SundanceUtils::Map<Set<int>, CellFilter> > 
      funcDomains(const Mesh& mesh,
                  const SundanceUtils::Map<CellFilter, Set<int> >& fmap,
                  SundanceUtils::Map<CellFilter, SundanceUtils::Map<Set<int>, CellSet> >& inputToChildrenMap);

      static SundanceUtils::Map<CellFilter, Set<int> > domainToFuncSetMap(const Array<Set<CellFilter> >& filters) ;

    private:

      static Set<CellFilter> reduceCellFilters(const Mesh& mesh,
                                               const Set<CellFilter>& inputSet) ;

      bool hasUnks() const ;

      bool unksAreOmnipresent() const ;

      bool testsAreOmnipresent() const ;

      bool regionIsMaximal(int r) const ;

      bool isSymmetric(int block) const ;

      void markBCRows(int block) ;

      const MPIComm& comm() const {return mesh().comm();}

      void init();

      Mesh mesh_;

      RefCountPtr<EquationSet> eqn_;

      Array<RefCountPtr<DOFMapBase> > rowMap_;

      Array<RefCountPtr<DOFMapBase> > colMap_;

      Array<RefCountPtr<Array<int> > > isBCRow_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
