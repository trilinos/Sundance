/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASICSIMPLICIALMESHTYPE_H
#define SUNDANCE_BASICSIMPLICIALMESHTYPE_H


#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceBasicSimplicialMesh.hpp"

namespace SundanceStdMesh
{
  using namespace Teuchos;
using namespace SundanceUtils;
  
  /**
   * BasicSimplicialMeshType is used to create
   * BasicSimplicialMesh objects.
   */
  class BasicSimplicialMeshType : public MeshTypeBase
  {
  public:
    /** Empty ctor */
    BasicSimplicialMeshType() {;}

    /** virtual dtor */
    virtual ~BasicSimplicialMeshType(){;}

    /** Create a mesh of the given dimension */
    virtual RefCountPtr<MeshBase> createEmptyMesh(int dim,
                                                  const MPIComm& comm) const 
    {return rcp(new BasicSimplicialMesh(dim, comm));}

    /** */
    string describe() const {return "BasicSimplicialMeshType";}

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshTypeBase> getRcp() {return rcp(this);}
#endif  /* DOXYGEN_DEVELOPER_ONLY */   
      
  };
}
#endif
