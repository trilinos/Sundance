#ifndef SUNDANCE_NULLCELLFILTER_BASEH
#define SUNDANCE_NULLCELLFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellFilterBase.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  
  namespace FrameworkInterface
  {
    /** 
     *
     * <h4> Notes for framework interface implementors </h4>
     *
     *  
     */
    class NullCellFilterBase : public CellFilterBase
    {
    public:
      /** Empty ctor */
      NullCellFilterBase();

      /** virtual dtor */
      virtual ~NullCellFilterBase(){;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const CellFilterBase* other) const ;

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return "NullCellFilterBase";}
      //@}

      /** */
      virtual RefCountPtr<CellFilterBase> getRcp() {return rcp(this);}

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
