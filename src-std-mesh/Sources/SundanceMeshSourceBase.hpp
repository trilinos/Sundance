/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHSOURCEBASE_H
#define SUNDANCE_MESHSOURCEBASE_H

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
     * MeshSourceBase provides the internal interface for mesh sources, i.e.,
     * objects such as mesh generators and file readers that can produce
     * a mesh object. The action of a mesh source should be independent
     * of the internal mesh representation used. To allow user-level
     * specification of the type of internal mesh representation to be
     * used, a MeshSourceBase is constructed with a MeshType object
     * which acts as a factory to produce empty meshes.
     *
     * Mesh sources are constructed with a MPI communicator. If the
     * communicator has more than one processor, the mesh created will
     * be distributed.
     *
     * <h4> Writing your own MeshSourceBase subtype </h4>
     *
     * The only method you will need to override is
     * <ul>
     * <li> <tt>virtual Mesh fillMesh() const </tt> 
     * </ul>
     * which is where you fill in a mesh object and return it. This method
     * should usually physically create the mesh with a call to createMesh(),
     * ensuring that the correct mesh representation type is created
     * using the MeshType factory with which the source was constructed.
     *
     * See the PartitionedLineMesher source code for a very simple
     * example of how to write a mesh source subtype. 
     *
     * Optionally, you can override the describe() method to 
     * provide more informative descriptive output than the string
     * <tt>"MeshSourceBase[unknown subtype]".</tt>
     */
    class MeshSourceBase : public TSFExtended::Handleable<MeshSourceBase>,
                           public TSFExtended::Printable,
                           public TSFExtended::Describable,
                           public Noncopyable,
                           public ObjectWithVerbosity<MeshSourceBase>
    {
    public:
      /** Construct with a mesh type and MPI communicator */
      MeshSourceBase(const MeshType& meshType,
                     const MPIComm& comm);

      /** virtual dtor */
      virtual ~MeshSourceBase(){;}

      
      /** Get a mesh from the source. If a mesh has already been created,
       * this method will return the cached mesh, otherwise it will 
       * build on with a call to fillMesh() */
      Mesh getMesh() const ;

      /** \name Extraction of attributes */
      //@{
      void getAttributes(RefCountPtr<Array<Array<double> > >& nodeAttributes,
                         RefCountPtr<Array<Array<double> > >& elemAttributes) const ;
      //@}

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << describe();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return "MeshSourceBase[unknown subtype]";}

      /** access to the MPI communicator */
      const MPIComm& comm() const {return comm_;}
      //@}
      
    protected:


      /** Get processor rank */
      int myRank() const {return comm().getRank();}
      /** Get number of processors */
      int nProc() const {return comm().getNProc();}

      /** Fill a mesh object with data from the source. Subclass
       * implementors will need to provide a fillMesh() for their
       * subclass. Implementors should use the createMesh() method
       * for the allocation of the new mesh. */
      virtual Mesh fillMesh() const = 0 ;

      /** createMesh() physically allocates the mesh object */
      Mesh createMesh(int dim) const ;

      /** internal access to the node attributes */
      RefCountPtr<Array<Array<double> > >& nodeAttributes() const 
      {return nodeAttributes_;}

      /** internal access to the element attributes */
      RefCountPtr<Array<Array<double> > >& elemAttributes() const 
      {return nodeAttributes_;}



    private:

      /** */
      mutable Mesh cachedMesh_;

      /** */
      mutable bool hasCachedMesh_;

      /** */
      MeshType meshType_;

      /** */
      MPIComm comm_;

      /** */
      mutable RefCountPtr<Array<Array<double> > > elemAttributes_;

      /** */
      mutable RefCountPtr<Array<Array<double> > > nodeAttributes_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
