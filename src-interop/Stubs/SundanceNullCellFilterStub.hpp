#ifndef SUNDANCE_NULLCELLFILTER_STUB_H
#define SUNDANCE_NULLCELLFILTER_STUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellFilterStub.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  
  namespace Internal
  {
    /** 
     *
     * <h4> Notes for framework interface implementors </h4>
     *
     *  
     */
    class NullCellFilterStub : public CellFilterStub
    {
    public:
      /** Empty ctor */
      NullCellFilterStub();

      /** virtual dtor */
      virtual ~NullCellFilterStub(){;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const CellFilterStub* other) const ;

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(std::ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string description() const 
      {return "NullCellFilterStub";}
      //@}

      /** */
      virtual RefCountPtr<CellFilterStub> getRcp() {return rcp(this);}

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
