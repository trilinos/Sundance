/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_HOMOGENEOUSDOFMAP_H
#define SUNDANCE_HOMOGENEOUSDOFMAP_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceDOFMapBase.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * A HomogeneousDOFMap is a DOF map for the special (and common)
     * case in which every function has the same basis and is defined
     * on every cell in the mesh. 
     */
    class HomogeneousDOFMap : public DOFMapBase
    {
    public:
      /** */
      HomogeneousDOFMap(const Mesh& mesh, 
                        const BasisFamily& basis,
                        int numFuncs);
                        


      /** */
      virtual ~HomogeneousDOFMap(){;}


      /** */
      virtual void getDOFsForCell(int cellDim, int cellLID,
                                  int funcID,
                                  Array<int>& dofs) const ;

      /** */
      virtual void print(ostream& os) const ;

    private:
      
      /** */
      bool hasBeenAssigned(int cellDim, int cellLID) const 
      {return dofs_[cellDim][cellLID][0] != uninitializedVal();}

      /** */
      void initMap();

      /** */
      void setDOFs(int cellDim, int cellLID, int& nextDOF);

      /** */
      void shareDOFs(int cellDim,
                     const Array<Array<int> >& outgoingCellRequests);

      /** */
      void computeOffsets(int dim, int localCount);

      static int uninitializedVal() {return -1;}

      int dim_;

      CellSet cells_;

      Array<Array<Array<int> > > dofs_;

      Array<int> funcID_;

      Array<Array<Array<Array<int> > > > localNodePtrs_;

      Array<int> nNodesPerCell_;

      Array<int> totalNNodesPerCell_;

      Array<Array<int> > numFacets_;

      bool basisIsContinuous_;
    };
  }
}

#endif
