/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHFILTER_H
#define SUNDANCE_MESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdMesh
{
  /**
   * MeshTransformation is the user-level interface for mesh filters, i.e.,
   * objects that take an input mesh and produce a new mesh. Examples
   * of filter operations are refinement, load balancing,
   * and extrusion from 2D to 3D. 
   *
   * <h4> Example: </h4> extrude a 2D mesh into 2D
   * \code
   * // create a 2D mesh 
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 10, 1,
   *                                                    0.0, 1.0, 10, 1,
   *                                                    meshType);
   * Mesh mesh2D = mesher.getMesh();
   * // create a filter for extruding 2 levels between z=0.0 and z=0.2
   * MeshTransformation extruder = new ExtrusionMeshTransformation(0.0, 0.2, 2);
   * // perform the extrusion
   * Mesh mesh3D = extruder.apply(mesh2D);
   * \endcode
   */
  class MeshTransformation : public TSFExtended::Handle<MeshTransformationBase>
  {
  public:
    /** Construct an empty mesh filter object */
    MeshTransformation();

    /** Construct from a raw pointer to a mesh filter subtype */
    MeshTransformation(Handleable<MeshTransformationBase>* rawPtr);

    /** Construct from a smart pointer to a mesh filter subtype */
    MeshTransformation(const RefCountPtr<MeshTransformationBase>& smartPtr);

    /** apply the filter to create a new mesh */
    Mesh apply(const Mesh& inputMesh) const ;

    const bool& serializeLocal() const {return serializeLocal_;}

    bool& serializeLocal() {return serializeLocal_;}
  private:
    bool serializeLocal_;
    
  };
}

#endif
