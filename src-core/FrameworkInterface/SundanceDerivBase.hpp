/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DERIVBASE_H
#define SUNDANCE_DERIVBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceExceptions.hpp"
#include "TSFHandleable.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY



namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      using namespace Teuchos;
      class Deriv;

      /**
       * DerivBase is a base class for first-order functional and spatial
       * differentiation in the region of internal AD. Multiple derivatives
       * are represented with the DerivSet object.
       */
      class DerivBase : public TSFExtended::Handleable<DerivBase>
        {
        public:
          /** */
          DerivBase();

          /** */
          virtual ~DerivBase(){;}

          /** */
          virtual bool lessThan(const Deriv& other) const = 0 ;

          /** */
          virtual string toString() const = 0 ;

        private:

        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */


#endif
