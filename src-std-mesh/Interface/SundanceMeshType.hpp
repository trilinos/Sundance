/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHTYPE_H
#define SUNDANCE_MESHTYPE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTypeBase.hpp"
#include "SundanceMesh.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdMesh
{
  using namespace Teuchos;
  using namespace SundanceUtils;
  using namespace TSFExtended;

  /**
   * Class MeshType is a user-level object for specification of which
   * internal mesh representation is to be used when building or reading
   * a mesh. An example of using a MeshType to control the creation 
   * of a mesh with a TriangleMeshReader is as follows: 
   * \code
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
   * \endcode
   * The internal representation of the mesh will be as a BasicSimplicialMesh
   * object. 
   */
  class MeshType : public TSFExtended::Handle<Internal::MeshTypeBase>
  {
  public:
    /** Construct an empty mesh type object */
    MeshType();

    /** Construct from a raw pointer to a mesh type subtype */
    MeshType(Handleable<Internal::MeshTypeBase>* rawPtr);

    /** Construct from a smart pointer to a mesh type subtype */
    MeshType(const RefCountPtr<Internal::MeshTypeBase>& smartPtr);

    /** Create a mesh of the given dimension */
    Mesh createEmptyMesh(int dim, const MPIComm& comm) const ;
    
  };
}

#endif
