/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CREATABLEMESH_H
#define SUNDANCE_CREATABLEMESH_H


#ifndef DOXYGEN_DEVELOPER_ONLY



#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     *
     */
    class CreatableMesh : public MeshBase
    {
    public:
      /** */
      CreatableMesh(int dim, const MPIComm& comm)
        : MeshBase(dim, comm) {;}

      /** */
      virtual ~CreatableMesh(){;}

      /** */
      virtual void estimateNumVertices(int nPts) {;}

      /** */
      virtual void estimateNumElements(int nElems) {;}

      

      /** */
      virtual int addVertex(int globalIndex, const Point& x,
                            int ownerProcID, int label) = 0 ;

      /** */
      virtual int addElement(int globalIndex, const Array<int>& vertLID,
                             int ownerProcID, int label) = 0 ;
  

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */


#endif
