/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHSOURCE_H
#define SUNDANCE_MESHSOURCE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdMesh
{
  /**
   * MeshSource is the user-level interface for objects such as
   * mesh generators and mesh file readers. A MeshSource can
   * create a mesh object with the <tt> getMesh() </tt> method,
   * and if node and element attributes are available, it can
   * access them with the getAttributes() method.
   *
   * Example: read input from a file celled "meshFile" 
   * in Shewchuk's Triangle format, and
   * create a mesh of type BasicSimplicialMesh distributed over the
   * MPI communicator MPI_COMM_WORLD.
   * \code
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
   * \endcode
   */
  class MeshSource : public TSFExtended::Handle<MeshSourceBase>
  {
  public:
    /** Construct an empty mesh source object */
    MeshSource();

    /** Construct from a raw pointer to a mesh source subtype */
    MeshSource(Handleable<MeshSourceBase>* rawPtr);

    /** Construct from a smart pointer to a mesh source subtype */
    MeshSource(const RefCountPtr<MeshSourceBase>& smartPtr);

    /** Create and return a mesh */
    Mesh getMesh() const ;

    /** Get any attributes associated with the nodes and elements in the
     * mesh. If no attributes exist, the arrays are empty. If the mesh
     * does not exist, it will be created with a cell to getMesh(). */
    void getAttributes(RefCountPtr<Array<Array<double> > >& nodeAttributes,
                       RefCountPtr<Array<Array<double> > >& elemAttributes) const ;

    static bool& serializeLocal() {static bool rtn=false; return rtn;}
  private:
    
  };
}

#endif
