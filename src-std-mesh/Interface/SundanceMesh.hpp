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

#ifndef SUNDANCE_MESH_H
#define SUNDANCE_MESH_H

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "SundanceIdentityReorderer.hpp"
#include "SundanceCellReorderer.hpp"
#include "SundanceHandle.hpp"

namespace SundanceStdMesh
{

using namespace Teuchos;
using namespace SundanceUtils;
using namespace Internal;

  /**
   * Mesh is the user-level object representing discrete geometry. 
   * The Mesh class is a handle to a MeshBase, which is an abstract interface
   * for meshes. 
   */
  class Mesh : public SundanceUtils::Handle<MeshBase>
  {
  public:

    /* */
    HANDLE_CTORS(Mesh, MeshBase);

    /** */
    int id() const {return ptr()->id();}

    /** 
     * Get the spatial dimension of the mesh
     */
    int spatialDim() const {return ptr()->spatialDim();}

    /** 
     * Get the number of cells having dimension dim
     */
    int numCells(int dim) const {return ptr()->numCells(dim);}

    /** 
     * Return the position of the i-th node
     */
    Point nodePosition(int i) const {return ptr()->nodePosition(i);}

    /** 
     * Return a view of the i-th node's position
     */
    const double* nodePositionView(int i) const {return ptr()->nodePositionView(i);}

    /** 
     * Return the centroid position of the cellLID-th cell of dimension
     * cellDim.
     */
    Point centroid(int cellDim, int cellLID) const 
    {return ptr()->centroid(cellDim, cellLID);}


    /** 
     * Get the outward normals for the batch of cells of dimension
     * spatialDim()-1. If any cell in the batch is not on the boundary,
     * an exception is thrown. 
     *
     * \param cellLIDs [in] LIDs for the cells whose normals are to be
     * computed. 
     * \param outwardNormals [out] Outward normal unit vectors for each
     * cell in the batch.
     */
    void outwardNormals(
      const Array<int>& cellLIDs,
      Array<Point>& outwardNormals
      ) const 
    {ptr()->outwardNormals(cellLIDs, outwardNormals);}

    /** 
     * Get a basis for the tangent space for each spatialDim()-1-dimensional
     * cell in a batch. The basis will be a set of spatialDim()-1 
     * spatialDim()-dimensional vectors represented as Point objects. 
     * The definition of the tangent basis is arbitrary up to a sign
     * in 2D or a rotation in 3D. In 2D, we will follow the sign convention
     * that cross(normal, tangent) = 1. In 3D, there is no simple
     * convention so the basis is chosen arbitrarily. Thus, in 3D it
     * is not a good idea to use the tangent basis in any way other than 
     * for specifying that both tangential components of some expression are
     * to be zero. 
     *
     * \param cellLIDs [in] LIDs for the cells whose tangent bases are to be
     * computed. 
     * \param tangentBasisVectors [out] Unit vectors spanning the tangent
     * space for each cell in the batch.
     */
    void tangentBasis(
      const Array<int>& cellLIDs,
      Array<Point>& tangentBasisVectors
      ) const ;
      
      
      

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
    void getJacobians(int cellDim, const Array<int>& cellLID,
                      CellJacobianBatch& jBatch) const 
    {ptr()->getJacobians(cellDim, cellLID, jBatch);}

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
                                  Array<double>& diameters) const 
    {ptr()->getCellDiameters(cellDim, cellLID, diameters);}


    /**
     * Map reference quadrature points to physical points on the
     * given cells. 
     */
    void pushForward(int cellDim, const Array<int>& cellLID,
                     const Array<Point>& refQuadPts,
                     Array<Point>& physQuadPts) const 
    {ptr()->pushForward(cellDim, cellLID, refQuadPts, physQuadPts);}

      

    /** 
     * Return the rank of the processor that owns the given cell
     */
    int ownerProcID(int cellDim, int cellLID) const 
    {return ptr()->ownerProcID(cellDim, cellLID);}
    

    /** 
     * Return the number of facets of the given cell
     */
    int numFacets(int cellDim, int cellLID, 
                  int facetDim) const 
    {return ptr()->numFacets(cellDim, cellLID, facetDim);}

    /** 
     * Return the local ID of a facet cell
     * @param cellDim dimension of the cell whose facets are being obtained
     * @param cellLID local index of the cell whose
     * facets are being obtained
     * @param facetDim dimension of the desired facet
     * @param facetIndex index into the list of the cell's facets
     */
    int facetLID(int cellDim, int cellLID,
                 int facetDim, int facetIndex,
                 int& facetOrientation) const 
    {return ptr()->facetLID(cellDim, cellLID, 
                            facetDim, facetIndex,
                            facetOrientation);}

    /**
     * Return by reference argument an
     *  array containing the LIDs of the facetDim-dimensional 
     * facets of the given cell
     */
    void getFacetArray(int cellDim, int cellLID, int facetDim, 
                       Array<int>& facetLIDs,
                       Array<int>& facetOrientations) const 
    {ptr()->getFacetArray(cellDim, cellLID, 
                          facetDim, facetLIDs,
                          facetOrientations);}

    /** 
     * Return a view of an element's zero-dimensional facets
     */
    const int* elemZeroFacetView(int cellLID) const 
    {return ptr()->elemZeroFacetView(cellLID);}

    /** 
     * Return by reference argument an array containing
     * the LIDs of the facetDim-dimensional facets of the
     * given batch of cells 
     */
    void getFacetLIDs(int cellDim, 
                      const Array<int>& cellLID,
                      int facetDim,
                      Array<int>& facetLID,
                      Array<int>& facetOrientations) const 
    {ptr()->getFacetLIDs(cellDim, cellLID, 
                         facetDim, facetLID, facetOrientations);}


    /** 
     * Return the number of maximal cofacets of the given cell
     */
    int numMaxCofacets(int cellDim, int cellLID) const 
    {return ptr()->numMaxCofacets(cellDim, cellLID);}

    /** 
     * Return the local ID of a maximal cofacet of a cell
     * @param cellDim dimension of the cell whose cofacets are being obtained
     * @param cellLID local index of the cell whose
     * cofacets are being obtained
     * @param cofacetIndex [in] index of the maximal cell 
     * into the list of the cell's cofacets
     * @param facetIndex [out] index of the calling cell
     * into the list of the maximal cell's facets
     */
    int maxCofacetLID(int cellDim, int cellLID,
                   int cofacetIndex,
                   int& facetIndex) const 
    {return ptr()->maxCofacetLID(cellDim, cellLID, cofacetIndex, facetIndex);}

    /** 
     * Get the LIDs of the maximal cofacets for a batch of 
     * codimension-one cells. 
     *
     * \param cellLIDs [in] array of LIDs of the cells whose cofacets are 
     * being obtained
     * \param cofacets [out] the batch of cofacets
     */
    void getMaxCofacetLIDs(const Array<int>& cellLIDs,
      MaximalCofacetBatch& cofacets) const 
      {ptr()->getMaxCofacetLIDs(cellLIDs, cofacets);}



    /** 
     * Find the cofacets of the given cell
     * @param cellDim dimension of the cell whose cofacets are being obtained
     * @param cellLID local index of the cell whose
     * cofacets are being obtained
     * @param cofacetDim dimension of the cofacets to get
     * @param cofacetLIDs LIDs of the cofacet
     */
    void getCofacets(int cellDim, int cellLID,
                     int cofacetDim, Array<int>& cofacetLIDs) const 
    {ptr()->getCofacets(cellDim, cellLID, cofacetDim, cofacetLIDs);}

    /** 
     * Find the local ID of a cell given its global index
     */
    int mapGIDToLID(int cellDim, int globalIndex) const 
    {return ptr()->mapGIDToLID(cellDim, globalIndex);}

    /** 
     * Determine whether a given cell GID exists on this processor
     */
    bool hasGID(int cellDim, int globalIndex) const 
    {return ptr()->hasGID(cellDim, globalIndex);}

    

    /** 
     * Find the global ID of a cell given its local index
     */
    int mapLIDToGID(int cellDim, int localIndex) const 
    {return ptr()->mapLIDToGID(cellDim, localIndex);}

    /**
     * Get the type of the given cell
     */
    CellType cellType(int cellDim) const 
    {return ptr()->cellType(cellDim);}

    /** Get the label of the given cell */
    int label(int cellDim, int cellLID) const 
    {return ptr()->label(cellDim, cellLID);}

    /** Get the labels for a batch of cells */
    void getLabels(int cellDim, const Array<int>& cellLID, Array<int>& labels) const 
    {ptr()->getLabels(cellDim, cellLID, labels);}

    /** Set the label for the given cell */
    void setLabel(int cellDim, int cellLID, int label)
    {ptr()->setLabel(cellDim, cellLID, label);}

    /** Get the list of all labels defined for cells of the given dimension */
    Set<int> getAllLabelsForDimension(int cellDim) const 
      {return ptr()->getAllLabelsForDimension(cellDim);}

    /** 
     * Return the number of labels associated with the given dimension.
     */
    virtual int numLabels(int cellDim) const 
      {return getAllLabelsForDimension(cellDim).size();}

    /** 
     * Get the cells associated with a specified label. The array 
     * cellLID will be filled with those cells of dimension cellDim
     * having the given label.
     */
    void getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const 
      {ptr()->getLIDsForLabel(cellDim, label, cellLIDs);}

    /** Get the MPI communicator over which the mesh is distributed */
    const MPIComm& comm() const {return ptr()->comm();}

    /** \name Mesh creation methods */
    //@{
    /** Allocate space for an estimated number of vertices */
    void estimateNumVertices(int nPts) 
    {creatableMesh()->estimateNumVertices(nPts);}
    
    /** Allocate space for an estimated number of elements */
    void estimateNumElements(int nElems) 
    {creatableMesh()->estimateNumElements(nElems);}

    
    /** Add a new vertex to the mesh */
    int addVertex(int globalIndex, const Point& x,
                  int ownerProcID, int label)
    {return creatableMesh()->addVertex(globalIndex, x, ownerProcID, label);}

    /** Add a new vertex to the mesh */
    int addElement(int globalIndex, const Array<int>& vertLID,
                  int ownerProcID, int label)
    {return creatableMesh()->addElement(globalIndex, vertLID, 
                                        ownerProcID, label);}

    
    /** \name Reordering */
    //@{
    /** Set the reordering method to be used with this mesh */
    void setReorderer(const CellReorderer& reorderer) 
    {ptr()->setReorderer(reorderer);}

    /** Set the reordering method to be used with this mesh */
    const CellReordererImplemBase* reorderer() const 
    {return ptr()->reorderer();}
    //@}
    
    /** */
    static CellReorderer& defaultReorderer()
    {
      static CellReorderer rtn = new IdentityReorderer();
      return rtn;
    }

    /** Work out global numberings for the cells of dimension cellDim */
    void assignIntermediateCellGIDs(int cellDim) 
    {
      if (!hasIntermediateGIDs(cellDim))
        ptr()->assignIntermediateCellGIDs(cellDim);
    }

    /** */
    bool hasIntermediateGIDs(int cellDim) const 
    {
      return ptr()->hasIntermediateGIDs(cellDim);
    }

    /** */
    void dump(const string& filename) const ;

    /** Test the consistency of the mesh numbering scheme 
     * across processors. This is meant as a check on Sundance's internal
     * logic rather than as a check on the validity of a user's mesh. */
    bool checkConsistency(const string& filename) const ;

    /** Test the consistency of the mesh numbering scheme 
     * across processors. This is meant as a check on Sundance's internal
     * logic rather than as a check on the validity of a user's mesh. */
    bool checkConsistency(ostream& os) const ;

  private:
    /** */
    IncrementallyCreatableMesh* creatableMesh();


    /** */
    bool checkVertexConsistency(ostream& os) const ;
    /** */
    bool checkCellConsistency(ostream& os, int dim) const ;

    /** */
    bool checkRemoteEntity(ostream& os, int p, int dim, int gid, 
                           int owner, bool mustExist, int& lid) const ;

    /** */
    bool testIdentity(ostream& os, int a, int b, const string& msg) const ;

    /** */
    bool testIdentity(ostream& os, 
                      const Array<int>& a,
                      const Array<int>& b, const string& msg) const ;
  };
}

#endif
