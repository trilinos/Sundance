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

#ifndef SUNDANCE_DOFMAPBASE_H
#define SUNDANCE_DOFMAPBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceBasisFamily.hpp"
#include "TSFObjectWithVerbosity.hpp"

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
     * 
     */
    class DOFMapBase : public TSFExtended::ObjectWithVerbosity<DOFMapBase>,
                       public TSFExtended::Printable
    {
    public:
      /** */
      DOFMapBase(const Mesh& mesh);

      
      /** */
      virtual ~DOFMapBase(){;}

      /** */
      bool isRemote(int cellDim, int cellLID, int& owner) const 
      {return (owner=mesh_.ownerProcID(cellDim, cellLID)) != localProcID_;}

      /** */
      virtual void getDOFsForCell(int cellDim, int cellLID,
                                  int funcID,
                                  Array<int>& dofs) const ;

      /** 
       * 
       */
      virtual void getDOFsForCellBatch(int cellDim, 
                                       const Array<int>& cellLID,
                                       Array<Array<int> >& dofs,
                                       Array<int>& nNodes) const = 0 ;

      /**
       * Given an input cell set, partition it into batches of cells
       * each having a homogeneous set of functions.
       *
       * \param inputCells [in] set of cells to be partitioned into batches
       * \param cellBatches [out] array of batches of cell indices.
       * \param homogSubregionIndices [out] indices to the function sets 
       * corresponding to each cell batch.
       *
       * Example: suppose the map has functions {0,1,2} on cells
       * 0 through 9, and functions {0,3} on cells 10 through 19.
       * Suppose the map stores function set {0,3} at index 0 and
       * function set {0,1,2} at index 1.
       * Calling this method with cellLID=[1,3,11,7,14] will result
       * in 
       * cellBatches = [[11,14], [1,3,7]]
       * homogSubregionIndices = [0,1]
       */
      virtual void 
      getHomogeneousCellBatches(const RefCountPtr<Array<int> >& inputCells,
                                Array<RefCountPtr<Array<int> > >& cellBatches,
                                Array<int>& homogSubregionIndices) const = 0 ;

      /** */
      virtual int numHomogeneousSubregions() const = 0 ;

      /** */
      virtual const Set<int>& getFuncSet(int homogSubregionIndex) const = 0 ;

      /** */
      virtual int chunkForFuncID(int homogSubregionIndex, int funcID) const = 0 ;

      /** */
      virtual int indexForFuncID(int homogSubregionIndex, int basisChunk, 
                                 int funcID) const = 0 ;

      /** */
      virtual int nFuncs(int homogSubregionIndex, int basisChunk) const = 0 ;

      /** */
      virtual int nBasisChunks(int homogSubregionIndex) const = 0 ;

      /** */
      virtual const BasisFamily& basis(int homogSubregionIndex, 
                                       int basisChunk) const = 0 ;

      /** */
      virtual const Array<int>& funcID(int homogSubregionIndex,
                                       int basisChunk) const = 0 ;

      /** */
      int lowestLocalDOF() const {return lowestLocalDOF_;}

      /** */
      bool isLocalDOF(int dof) const 
      {return (dof >= lowestLocalDOF_ && dof < lowestLocalDOF_+numLocalDOFs_);}

      /** */
      int numLocalDOFs() const {return numLocalDOFs_;}

      /** */
      int numDOFs() const {return numDOFs_;}

      /** */
      const RefCountPtr<Array<int> >& ghostIndices() const 
      {return ghostIndices_;}

      /** */
      virtual void print(ostream& os) const = 0 ;


    protected:
      
      void verifySubregionIndex(int subregionIndex) const ;

      void verifyBasisChunkIndex(int subregionIndex,
                                 int basisChunk) const ;

      void setLowestLocalDOF(int low) {lowestLocalDOF_ = low;}

      void setNumLocalDOFs(int numDOFs) {numLocalDOFs_ = numDOFs;}

      void setTotalNumDOFs(int numDOFs) {numDOFs_ = numDOFs;}

      const Mesh& mesh() const {return mesh_;}

      const MPIComm& comm() const {return mesh().comm();}

      void addGhostIndex(int dof) {ghostIndices_->append(dof);}

      static Time& dofLookupTimer() ;

      static Time& batchedDofLookupTimer() ;

    private:
      int localProcID_;

      Mesh mesh_;

      int lowestLocalDOF_;

      int numLocalDOFs_;

      int numDOFs_;

      RefCountPtr<Array<int> > ghostIndices_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
