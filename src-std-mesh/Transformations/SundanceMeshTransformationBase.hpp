/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHFILTERBASE_H
#define SUNDANCE_MESHFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "TSFHandleable.hpp"
#include "TSFDescribable.hpp"
#include "TSFPrintable.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceCreatableMesh.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     * MeshSourceBase provides the internal interface for mesh filters, i.e.,
     * objects that take an input mesh and produce a new mesh. Examples
     * of filter operations are refinement, load balancing,
     * and extrusion from 2D to 3D. 
     * The action of a mesh filter should be independent
     * of the internal mesh representation used. To allow user-level
     * specification of the type of internal mesh representation to be
     * used, a MeshTransformationBase is constructed with a MeshType object
     * which acts as a factory to produce empty output meshes.
     *
     * If the
     * communicator has more than one processor, the mesh created will
     * be distributed.
     *
     * <h4> Writing your own MeshTransformationBase subtype </h4>
     *
     * The only method you will need to override is
     * <ul>
     * <li> <tt>virtual Mesh apply(const Mesh& inputMesh) const = 0 </tt> 
     * </ul>
     * which is where you do the filter action and return an output
     * mesh. This method
     * should usually physically create the mesh with a call to createMesh(),
     * ensuring that the correct mesh representation type is created
     * using the MeshType factory with which the filter was constructed.
     *
     * See the ExtrustionMeshTransformation source code for a very simple
     * example of how to write a mesh filter subtype. 
     *
     * Optionally, you can override the describe() method to 
     * provide more informative descriptive output than the string
     * <tt>"MeshTransformationBase[unknown subtype]".</tt>
     */
    class MeshTransformationBase : public TSFExtended::Handleable<MeshTransformationBase>,
                           public TSFExtended::Printable,
                           public TSFExtended::Describable,
                           public Noncopyable,
                           public TSFExtended::ObjectWithVerbosity<MeshTransformationBase>
    {
    public:
      /** Construct with a mesh type, which specifies
       *  the type of mesh to be built when the filter is applied. */
      MeshTransformationBase(const MeshType& meshType)
        : meshType_(meshType) {;}

      /** virtual dtor */
      virtual ~MeshTransformationBase(){;}

      
      /** Apply the filter to the given input mesh, 
       *  producing an output mesh */
      virtual Mesh apply(const Mesh& inputMesh) const = 0 ;

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << describe();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return "MeshTransformationBase[unknown subtype]";}
      //@}
      
    protected:

      /** createMesh() allocates the mesh object with a call to 
       * meshType's createMesh() method. */
      Mesh createMesh(int dim, const MPIComm& comm) const ;

    private:
      /** */
      MeshType meshType_;


    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif