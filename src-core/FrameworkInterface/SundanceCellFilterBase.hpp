/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLFILTERBASE_H
#define SUNDANCE_CELLFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceObjectWithInstanceID.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "Teuchos_XMLObject.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  
  namespace FrameworkInterface
  {
    /** 
     * Base class for cell filter objects, i.e., objects that can 
     * select a subset of mesh cells on which an integral or 
     * BC is to be applied.
     *
     * <h4> Notes for framework interface implementors </h4>
     *
     *  
     */
    class CellFilterBase : public TSFExtended::Handleable<CellFilterBase>,
                           public TSFExtended::Printable,
                           public TSFExtended::Describable,
                           public Noncopyable,
                           public ObjectWithInstanceID<CellFilterBase>
    {
    public:
      /** Empty ctor */
      CellFilterBase();

      /** virtual dtor */
      virtual ~CellFilterBase(){;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const CellFilterBase* other) const 
      {return id() < other->id();}

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return "CellFilterBase[id=" + Teuchos::toString(id()) + "]";}
      //@}

      /** */
      virtual RefCountPtr<CellFilterBase> makeNullDomain() const ;

      /** */
      virtual RefCountPtr<CellFilterBase> getRcp() {return rcp(this);}

      /** */
      bool isNullDomain() const ;

      

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
