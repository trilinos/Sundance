/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHCREATIONINTERFACE_H
#define SUNDANCE_MESHCREATIONINTERFACE_H


#ifndef DOXYGEN_DEVELOPER_ONLY



#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
  namespace Internal
  {
    /**
     *
     */
    class MeshCreationInterface : public MeshBase
    {
    public:
      /** */
      MeshCreationInterface(int dim, const MPIComm& comm)
        : MeshBase(dim, comm) {;}

      /** */
      virtual ~MeshCreationInterface(){;}

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
