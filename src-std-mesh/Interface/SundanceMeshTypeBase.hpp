/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHTYPEBASE_H
#define SUNDANCE_MESHTYPEBASE_H


#include "SundanceDefs.hpp"
#include "TSFHandleable.hpp"
#include "TSFDescribable.hpp"
#include "TSFPrintable.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceMeshBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  namespace Internal
  {
    using namespace Teuchos;

    /**
     * MeshTypeBase is a factory class for empty meshes, allowing generic
     * mesh sources to build a mesh of user-specified type. 
     */
    class MeshTypeBase : public TSFExtended::Handleable<MeshTypeBase>,
                         public TSFExtended::Printable,
                         public TSFExtended::Describable,
                         public Noncopyable
    {
    public:
      /** Empty ctor */
      MeshTypeBase() {;}

      /** virtual dtor */
      virtual ~MeshTypeBase(){;}

      /** Create a mesh of the given dimension */
      virtual RefCountPtr<MeshBase> createEmptyMesh(int dim,
                                                    const MPIComm& comm) const = 0 ;

      /** */
      virtual void print(ostream& os) const {os << describe();}
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */   

#endif
