/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ZEROEXPR_H
#define SUNDANCE_ZEROEXPR_H

#include "SundanceConstantExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      class ZeroExpr : public ConstantExpr
        {
        public:
          /** */
          ZeroExpr();

          /** */
          virtual ~ZeroExpr() {;}

          /** */
          virtual bool hasNonzeroDeriv(const MultipleDeriv& /* d */) const 
          {return false;}
          
          /** */
          virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
        protected:
        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
