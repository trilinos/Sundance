/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLREORDERERBASE_H
#define SUNDANCE_CELLREORDERERBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceCellReordererImplemBase.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    using namespace Teuchos;
using namespace SundanceUtils;

    /**
     * Factory class to instantiate cell reorderers for specific meshes
     */
    class CellReordererFactoryBase 
      : public TSFExtended::Handleable<CellReordererFactoryBase>, 
        public Noncopyable,
        public TSFExtended::Printable,
        public TSFExtended::Describable
    {
    public:
      /** */
      CellReordererFactoryBase() {;}
      
      /** virtual dtor */
      virtual ~CellReordererFactoryBase(){;}

      /** */
      virtual string describe() const {return typeid(*this).name();}

      /** */
      virtual void print(ostream& os) const {os << describe();}

      /** Instantiate a factory */
      virtual RefCountPtr<CellReordererImplemBase> 
      createInstance(const MeshBase* mesh) const = 0 ;
    };

    /**
     * Factory class to instantiate cell reorderers for specific meshes
     */
    template <class T>
    class GenericCellReordererFactory : public CellReordererFactoryBase
    {
    public:
      /** */
      GenericCellReordererFactory() {;}
      
      /** virtual dtor */
      virtual ~GenericCellReordererFactory(){;}

      /** Instantiate a factory */
      virtual RefCountPtr<CellReordererImplemBase> 
      createInstance(const MeshBase* mesh) const 
      {
        return rcp(new T(mesh));
      }
      
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
