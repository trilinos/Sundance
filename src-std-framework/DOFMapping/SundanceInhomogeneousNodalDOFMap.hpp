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

#ifndef SUNDANCE_INHOMOGENEOUSNODALDOFMAP_H
#define SUNDANCE_INHOMOGENEOUSNODALDOFMAP_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using SundanceUtils::Map;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * 
     */
    class InhomogeneousNodalDOFMap : public DOFMapBase
    {
    public:
      /** */
      InhomogeneousNodalDOFMap(const Mesh& mesh, 
        const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap, const ParameterList& verbParams = *DOFMapBase::defaultVerbParams());
      
      /** */
      virtual ~InhomogeneousNodalDOFMap(){;}

      /** */
      RefCountPtr<const MapStructure> 
      getDOFsForCellBatch(int cellDim,
                          const Array<int>& cellLID,
                          const Set<int>& requestedFuncSet,
                          Array<Array<int> >& dofs,
                          Array<int>& nNodes) const ;

      /** */
      void getFunctionDofs(int cellDim,
                           const Array<int>& cellLID,
                           const Array<int>& facetLID,
                           const Array<int>& funcs,
                           Array<Array<int> >& dofs) const ;

      /** */
      RefCountPtr<const Set<int> >
      allowedFuncsOnCellBatch(int cellDim,
                              const Array<int>& cellLID) const ;

      /** */
      const Array<CellFilter>& funcDomains() const {return funcDomains_;}

      /** */
      virtual void print(ostream& os) const ;


    protected:

      /** */
      Array<int> dofsOnCell(int cellDim, int cellLID, const Set<int>& reqFuncs) const ;
                              

      void init();

      void computeOffsets(int localCount)  ;

      void shareRemoteDOFs(const Array<Array<int> >& remoteNodes);

      void assignNode(int fLID,
                      int funcComboIndex,
                      int dofOffset,
                      int nFuncs,
                      Array<Array<int> >& remoteNodes,
                      Array<int>& hasProcessedCell,
                      int& nextDOF) ;

      int dim_;
      RCP<BasisDOFTopologyBase> basis_;
      int nTotalFuncs_;
      Array<CellFilter> funcDomains_;

      Array<Array<int> > nodeDofs_;
      Array<Array<int> > elemDofs_;
      Array<int> nodeToFuncSetIndexMap_;
      Array<int> elemToFuncSetIndexMap_;
      Array<Set<int> > elemFuncSets_;
      Array<Set<int> > nodalFuncSets_;
      Array<int> nodeToOffsetMap_;
      Array<int> elemToOffsetMap_;

      Array<Array<int> > funcIndexWithinNodeFuncSet_;

      Array<RefCountPtr<const MapStructure> > elemStructure_;
      Array<RefCountPtr<const MapStructure> > nodeStructure_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
