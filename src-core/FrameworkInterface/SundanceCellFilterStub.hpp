/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLFILTERSTUB_H
#define SUNDANCE_CELLFILTERSTUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceObjectWithInstanceID.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "Teuchos_XMLObject.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  
  namespace Internal
  {
    /** 
     * Stub class for cell filter objects, i.e., objects that can 
     * select a subset of mesh cells on which an integral or 
     * BC is to be applied.
     *
     * <h4> Notes for framework interface implementors </h4>
     *
     *  
     */
    class CellFilterStub : public TSFExtended::Handleable<CellFilterStub>,
                           public TSFExtended::Printable,
                           public TSFExtended::Describable,
                           public Noncopyable,
                           public ObjectWithInstanceID<CellFilterStub>,
                           public TSFExtended::ObjectWithVerbosity<CellFilterStub>
    {
    public:
      /** Empty ctor */
      CellFilterStub();

      /** virtual dtor */
      virtual ~CellFilterStub(){;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const CellFilterStub* other) const 
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
      {return "CellFilterStub[id=" + Teuchos::toString(id()) + "]";}
      //@}

      /** */
      virtual RefCountPtr<CellFilterStub> makeNullRegion() const ;

      /* */
      GET_RCP(CellFilterStub);

      /** */
      bool isNullRegion() const ;

      

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif