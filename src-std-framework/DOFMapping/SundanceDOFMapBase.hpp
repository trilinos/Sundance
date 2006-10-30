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
#include "SundanceMapStructure.hpp"
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
      virtual RefCountPtr<const MapStructure> 
      getDOFsForCellBatch(int cellDim,
                          const Array<int>& cellLID,
                          const Set<int>& requestedFuncSet,
                          Array<Array<int> >& dofs,
                          Array<int>& nNodes) const = 0 ;

      /** */
      virtual RefCountPtr<const Set<int> >
      allowedFuncsOnCellBatch(int cellDim,
                              const Array<int>& cellLID) const = 0 ;

      /** */
      virtual const Array<CellFilter>& funcDomains() const = 0 ;

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

      /** */
      virtual bool isHomogeneous() const {return false;}


    protected:
      

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
