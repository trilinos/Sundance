/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DERIVSET_H
#define SUNDANCE_DERIVSET_H



#include "SundanceDefs.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceSet.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {

      /**
       * DerivSet is a set of multiple derivatives, sorted by 
       * increasing differentiation order.
       *
       * The sorting is important because we should evaluate derivatives
       * of products in decreasing order of differentiation order.
       */
      class DerivSet : public SundanceUtils::Set<MultipleDeriv,
        increasingOrder<MultipleDeriv> >
        {
        public:
          /** Construct an empty set of functional derivatives */
          DerivSet();

        };
    }
}

namespace Teuchos
{
  using namespace SundanceUtils;
  using namespace SundanceCore::Internal;
  using namespace SundanceCore;
  inline string toString(const SundanceCore::Internal::DerivSet& d)
  {
    return d.toString();
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
