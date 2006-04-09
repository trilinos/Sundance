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

#ifndef SUNDANCE_MESHBASE_H
#define SUNDANCE_MESHBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "SundanceCellType.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceObjectWithInstanceID.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceCellReorderer.hpp"
#include "SundanceCellReordererImplemBase.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    using namespace Teuchos;
using namespace SundanceUtils;
    class CellJacobianBatch;

    /**
     * Abstract interface for meshes.
     *
     */
    class MeshBase : public TSFExtended::ObjectWithVerbosity<MeshBase>,
                     public ObjectWithInstanceID<MeshBase>
    {
    public:
      /** */
      MeshBase(int dim, const MPIComm& comm);

      /** */
      virtual ~MeshBase(){;}

      /** 
       * Get the spatial dimension of the mesh
       */
      int spatialDim() const {return dim_;}

      /** 
       * Get the number of cells having dimension dim
       */
      virtual int numCells(int dim) const = 0 ;

      /** 
       * Return the position of the i-th node
       */
      virtual Point nodePosition(int i) const = 0 ;

      /** 
       * Compute the jacobians of a batch of cells, returning the 
       * result via reference argument
       *
       * @param cellDim dimension of the cells whose Jacobians are to
       * be computed
       * @param cellLID local indices of the cells for which Jacobians
       * are to be computed
       * @param jBatch reference to the resulting Jacobian batch
       */
      virtual void getJacobians(int cellDim, const Array<int>& cellLID,
                                CellJacobianBatch& jBatch) const { ; }

      /** 
       * Compute the diameters of a batch of cells,
       * result via reference argument
       *
       * @param cellDim dimension of the cells whose diameters are to
       * be computed
       * @param cellLID local indices of the cells for which diameters
       * are to be computed
       * @param diameters reference to the array of cell diameters
       */
      virtual void getCellDiameters(int cellDim, const Array<int>& cellLID,
                                    Array<double>& diameters) const { ; }

      /**
       * Map reference quadrature points to physical points on the
       * given cells. 
       */
      virtual void pushForward(int cellDim, const Array<int>& cellLID,
                               const Array<Point>& refQuadPts,
                               Array<Point>& physQuadPts) const { ; }

      

      /** 
       * Return the rank of the processor that owns the given cell
       */
      virtual int ownerProcID(int cellDim, int cellLID) const = 0 ;

      /** 
       * Return the number of facets of the given cell
       */
      virtual int numFacets(int cellDim, int cellLID, 
                            int facetDim) const = 0 ;
 
      /** 
       * Return the local ID of a facet cell
       * @param cellDim dimension of the cell whose facets are being obtained
       * @param cellLID local index of the cell whose
       * facets are being obtained
       * @param facetDim dimension of the desired facet
       * @param facetIndex index into the list of the cell's facets
       */
      virtual int facetLID(int cellDim, int cellLID,
                           int facetDim, int facetIndex,
                           int& facetOrientation) const = 0 ;

      /** 
       * Return by reference argument an array containing
       * the LIDs of the facetDim-dimensional facets of the
       * given batch of cells 
       */
      virtual void getFacetLIDs(int cellDim, 
                                const Array<int>& cellLID,
                                int facetDim,
                                Array<int>& facetLID,
                                Array<int>& facetOrientations) const = 0 ;

      /**
       * Return by reference argument an
       * array containing the LIDs of the facetDim-dimensional 
       * facets of the given cell. The default implementation 
       * loops over calls to facetLID(). Subclasses can 
       * provide a more efficient implementation if desired. 
       */
      void getFacetArray(int cellDim, int cellLID, int facetDim, 
                         Array<int>& facetLIDs,
                         Array<int>& facetOrientations) const ;


      /** 
       * Return the number of maximal cofacets of the given cell
       */
      virtual int numCofacets(int cellDim, int cellLID) const = 0 ;

      /** 
       * Return the local ID of a maximal cofacet of a cell
       * @param cellDim dimension of the cell whose cofacets are being obtained
       * @param cellLID local index of the cell whose
       * cofacets are being obtained
       * @param cofacetIndex index into the list of the cell's facets
       */
      virtual int cofacetLID(int cellDim, int cellLID,
                             int cofacetIndex, 
                             int& facetIndex) const = 0 ;

      /** 
       * Find the local ID of a cell given its global index
       */
      virtual int mapGIDToLID(int cellDim, int globalIndex) const = 0 ;

    
      /** 
       * Determine whether a given cell GID exists on this processor
       */
      virtual bool hasGID(int cellDim, int globalIndex) const = 0 ;


      /** 
       * Find the global ID of a cell given its local index
       */
      virtual int mapLIDToGID(int cellDim, int localIndex) const = 0 ;

      /**
       * Get the cell type used in the given dimension.
       */
      virtual CellType cellType(int cellDim) const = 0 ;


      /** Get the label of the given cell */
      virtual int label(int cellDim, int cellLID) const = 0 ;

      /** Set the label for the given cell */
      virtual void setLabel(int cellDim, int cellLID, int label) = 0 ;

      /** Work out global numberings for the cells of dimension cellDim */
      virtual void assignIntermediateCellGIDs(int cellDim) = 0 ;

      /** */
      virtual bool hasIntermediateGIDs(int dim) const = 0 ;


      /** Return the MPI communicator over which this mesh is
       * distributed */
      const MPIComm& comm() const {return comm_;}

      
      /** \name Reordering */
      //@{
      /** Set the reordering method to be used with this mesh */
      void setReorderer(const CellReorderer& reorderer) 
      {reorderer_ = reorderer.createInstance(this);}

      /** Set the reordering method to be used with this mesh */
      const CellReordererImplemBase* reorderer() const 
      {return reorderer_.get();}
      //@}


      /** Whether to stagger output in parallel. Set to true for
       * readable output in parallel debugging sessions. This should
       * be normally, as it causes one synchronization point per process. */
      static bool& staggerOutput() {static bool rtn=false; return rtn;}
    
    private:
      int dim_;

      MPIComm comm_;

      RefCountPtr<CellReordererImplemBase> reorderer_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
