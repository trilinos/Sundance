/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESH_H
#define SUNDANCE_MESH_H

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceCreatableMesh.hpp"
#include "SundanceIdentityReorderer.hpp"
#include "SundanceCellReorderer.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
using namespace SundanceUtils;
  using namespace Internal;
  /**
   *
   */
  class Mesh : public Handle<MeshBase>
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
                 int facetDim, int facetIndex) const 
    {return ptr()->facetLID(cellDim, cellLID, facetDim, facetIndex);}

    /**
     * Return by reference argument an
     *  array containing the LIDs of the facetDim-dimensional 
     * facets of the given cell
     */
    void getFacetArray(int cellDim, int cellLID, int facetDim, 
                       Array<int>& facetLIDs) const 
    {ptr()->getFacetArray(cellDim, cellLID, facetDim, facetLIDs);}

    /** 
     * Return by reference argument an array containing
     * the LIDs of the facetDim-dimensional facets of the
     * given batch of cells 
     */
    void getFacetLIDs(int cellDim, 
                      const Array<int>& cellLID,
                      int facetDim,
                      Array<int>& facetLID) const 
    {ptr()->getFacetLIDs(cellDim, cellLID, facetDim, facetLID);}


    /** 
     * Return the number of maximal cofacets of the given cell
     */
    int numCofacets(int cellDim, int cellLID) const 
    {return ptr()->numCofacets(cellDim, cellLID);}

    /** 
     * Return the local ID of a maximal cofacet of a cell
     * @param cellDim dimension of the cell whose cofacets are being obtained
     * @param cellLID local index of the cell whose
     * cofacets are being obtained
     * @param cofacetIndex index into the list of the cell's facets
     */
    int cofacetLID(int cellDim, int cellLID,
                   int cofacetIndex) const 
    {return ptr()->cofacetLID(cellDim, cellLID, cofacetIndex);}

    /** 
     * Find the local ID of a cell given its global index
     */
    int mapGIDToLID(int cellDim, int globalIndex) const 
    {return ptr()->mapGIDToLID(cellDim, globalIndex);}

    

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
    void assignIntermediateCellOwners(int cellDim) 
    {
      ptr()->assignIntermediateCellOwners(cellDim);
    }

    /** */
    bool hasIntermediateGIDs(int cellDim) const 
    {
      return ptr()->hasIntermediateGIDs(cellDim);
    }
    

  private:
    /** */
    CreatableMesh* creatableMesh();
  };
}

#endif
