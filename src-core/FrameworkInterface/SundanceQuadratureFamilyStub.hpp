/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREFAMILYSTUB_H
#define SUNDANCE_QUADRATUREFAMILYSTUB_H

#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_XMLObject.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    using namespace Teuchos;
    using namespace Internal;

    class QuadratureFamilyStub 
      : public TSFExtended::Handleable<QuadratureFamilyStub>,
        public TSFExtended::Printable,
        public TSFExtended::Describable,
        public Noncopyable
    {
    public:
      /** Empty ctor */
      QuadratureFamilyStub(int order);

      /** virtual dtor */
      virtual ~QuadratureFamilyStub(){;}

      /** Return the order of the quadrature rule */
      int order() const {return order_;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const QuadratureFamilyStub* other) const 
      {return order() < other->order();}

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return "QuadratureFamilyStub[order=" + Teuchos::toString(order()) 
         +  "]";}
      //@}

      /** */
      virtual RefCountPtr<QuadratureFamilyStub> getRcp() {return rcp(this);}

      /** */
      static RefCountPtr<QuadratureFamilyStub>& defaultQuadrature()
      {static RefCountPtr<QuadratureFamilyStub> rtn; return rtn;}
      
    private:
      int order_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
