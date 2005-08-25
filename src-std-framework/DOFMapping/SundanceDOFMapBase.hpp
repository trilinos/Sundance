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

      /** */
      virtual void getDOFsForCellBatch(int cellDim, const Array<int>& cellLID,
                                       Array<int>& dofs, 
                                       unsigned int& nNodes) const = 0 ;

      /** */
      const CellSet& cellSet(int i) const {return cellSets_[i];}

      /** */
      int numCellSets() const {return cellSets_.size();}

      /** */
      const Array<int>& funcIDOnCellSet(int i) const 
      {return funcIDOnCellSets_[i];}

      /** */
      int cellDimOnCellSet(int i) const 
      {return cellDimOnCellSets_[i];}

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


    protected:

      void setLowestLocalDOF(int low) {lowestLocalDOF_ = low;}

      void setNumLocalDOFs(int numDOFs) {numLocalDOFs_ = numDOFs;}

      void setTotalNumDOFs(int numDOFs) {numDOFs_ = numDOFs;}

      const Mesh& mesh() const {return mesh_;}

      const MPIComm& comm() const {return mesh().comm();}

      Array<CellSet>& cellSets() {return cellSets_;}

      Array<Array<int> >& funcIDOnCellSets() {return funcIDOnCellSets_;}

      Array<int>& cellDimOnCellSets() {return cellDimOnCellSets_;}

      void addGhostIndex(int dof) {ghostIndices_->append(dof);}

    private:
      int localProcID_;

      Mesh mesh_;

      Array<CellSet> cellSets_;

      Array<Array<int> > funcIDOnCellSets_;

      Array<int> cellDimOnCellSets_;

      int lowestLocalDOF_;

      int numLocalDOFs_;

      int numDOFs_;

      RefCountPtr<Array<int> > ghostIndices_;

      Array<Array<int> > dofsHaveBeenAssigned_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
