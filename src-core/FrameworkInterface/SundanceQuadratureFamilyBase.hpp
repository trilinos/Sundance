/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_QUADRATUREFAMILYBASE_H
#define SUNDANCE_QUADRATUREFAMILYBASE_H

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
  namespace FrameworkInterface
  {
    using namespace Teuchos;
    using namespace Internal;

    class QuadratureFamilyBase 
      : public TSFExtended::Handleable<QuadratureFamilyBase>,
        public TSFExtended::Printable,
        public TSFExtended::Describable,
        public Noncopyable
    {
    public:
      /** Empty ctor */
      QuadratureFamilyBase(int order);

      /** virtual dtor */
      virtual ~QuadratureFamilyBase(){;}

      /** Return the order of the quadrature rule */
      int order() const {return order_;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const QuadratureFamilyBase* other) const 
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
      {return "QuadratureFamilyBase[order=" + Teuchos::toString(order()) 
         +  "]";}
      //@}

      /** */
      virtual RefCountPtr<QuadratureFamilyBase> getRcp() {return rcp(this);}
      
    private:
      int order_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
