/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_COORDDERIV_H
#define SUNDANCE_COORDDERIV_H

#include "SundanceDefs.hpp"
#include "SundanceDerivBase.hpp"
#include "SundanceMultiIndex.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      using namespace Teuchos;

      /**
       * CoordDeriv represents a single first-order derivative with respect
       * to a coordinate function. It is used internally
       * to represent differentiation of
       * explicit coordinate dependence.
       *
       * @see Deriv
       *
       */
      class CoordDeriv : public DerivBase
        {
        public:
          /** */
          CoordDeriv(int dir);

          /** */
          virtual ~CoordDeriv(){;}

          /** */
          int dir() const {return dir_;}

          /** */
          virtual bool lessThan(const Deriv& other) const ;

          /** */
          virtual string toString() const ;

          /* handleable boilerplate */
          GET_RCP(DerivBase);

          /** */
          static int maxDim() {static int rtn=3; return rtn;}
        private:
          int dir_;
        };
    }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
