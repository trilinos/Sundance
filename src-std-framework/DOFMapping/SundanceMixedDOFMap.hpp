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

#ifndef SUNDANCE_MIXEDDOFMAP_H
#define SUNDANCE_MIXEDDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceDOFMapBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * A MixedDOFMap is a DOF map for the case where 
     * every function is defined
     * on every cell in the mesh, but where functions may have different bases. 
     */
    class MixedDOFMap : public DOFMapBase
    {
    public:
      /** */
      MixedDOFMap(const Mesh& mesh, 
                  const BasisArray& basis);
                        
      /** */
      virtual ~MixedDOFMap(){;}


     

      

      /** */
      virtual void getDOFsForCellBatch(int cellDim, 
                                       const Array<int>& cellLID,
                                       Array<Array<int> >& dofs,
                                       Array<int>& nNodes) const ;
      

      /** */
      virtual void print(ostream& os) const ;

    private:

      /** */
      inline int getInitialDOFForCell(int cellDim, int cellLID, int basisChunk) const
      {
        return dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]];
      }

      inline int* getInitialDOFPtrForCell(int cellDim, int cellLID, int basisChunk)
      {
        return &(dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]]);
      }

      inline const int* getInitialDOFPtrForCell(int cellDim, int cellLID, 
                                                int basisChunk) const 
      {
        return &(dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]]);
      }

      /** */
      void allocate(const Mesh& mesh);
      
      /** */
      void buildMaximalDofTable() const ;

      bool hasBeenAssigned(int cellDim, int cellLID) const 
      {return hasBeenAssigned_[cellDim][cellLID];}

      void markAsAssigned(int cellDim, int cellLID)
      {hasBeenAssigned_[cellDim][cellLID] = true;}

      /** */
      void initMap();

      /** */
      void setDOFs(int basisChunk, int cellDim, int cellLID, 
                   int& nextDOF, bool isRemote=false);

      /** */
      void shareDOFs(int cellDim,
                     const Array<Array<int> >& outgoingCellRequests);

      /** */
      void computeOffsets(int dim, int localCount);

      /** */
      const Array<int>& funcIDList() const {return funcIDOnCellSet(0);}

      /** */
      static int uninitializedVal() {return -1;}

      /** Table of bases for each chunk */
      Array<BasisFamily> chunkBasis_;

      /** Table of lists of funcID for each chunk */
      Array<Array<int> > chunkFuncIDs_;

      /** Map for lookup of chunk number given a function ID */
      Array<int> funcIDToChunkMap_;

      /** spatial dimension */
      int dim_;

      /** Tables of DOFs, indexed by dimension and chunk number.
       *
       * dof(cellDim, cellLID, chunk, func, node) 
       * = dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node]
       */
      Array<Array<Array<int> > > dofs_;

      /** DOFs for maximal cells, indexed by basis chunk number 
       *
       * dof(cellLID, chunk, func, node) 
       * = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node];
       */
      mutable Array<Array<int> > maximalDofs_;

      /** whether maximal DOFs have been tabulated */
      mutable bool haveMaximalDofs_;

      /** 
       * localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
       */
      Array<Array<Array<Array<Array<int> > > > > localNodePtrs_;

      /** The number of nodes per cell, for each basis function type, 
       * <i>not</i> including the nodes of the facets of the cell. Indexed as 
       * nNodesPerCell_[basis][dimension] */
      Array<Array<int> > nNodesPerCell_;

      /** The number of DOFs per cell, for each basis function type, 
       * <i>not</i> including the DOFs of the facets of the cell. Indexed as 
       * nDofsPerCell_[basis][dimension] */
      Array<Array<int> > nDofsPerCell_;

      /** The number of nodes per cell, for each basis function type, including
       * the nodes of the facets of the cell. Indexed as 
       * nNodesPerCell_[basis][dimension] */
      Array<Array<int> > totalNNodesPerCell_;

      /** The number of DOFs per cell, for each basis function type, including
       * the DOFs of the facets of the cell. Indexed as 
       * nDofsPerCell_[basis][dimension] */
      Array<Array<int> > totalNDofsPerCell_;

      /** Indicates whether the cells of each dimension have any DOFs in this
       * map, for any chunk. */
      Array<int> cellHasAnyDOFs_;

      /** number of facets of dimension facetDim for cells of dimension cellDim.
       * Indexed as numFacets_[cellDim][facetDim]
       */
      Array<Array<int> > numFacets_;

      /** Orientation of each edge or face as seen by the maximal cell
       * from which its DOFs were originally assigned. */
      Array<Array<int> > originalFacetOrientation_;

      /** */
      Array<Array<int> > hasBeenAssigned_;

    };
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  


#endif
