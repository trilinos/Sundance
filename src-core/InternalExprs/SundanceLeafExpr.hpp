/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LEAFEXPR_H
#define SUNDANCE_LEAFEXPR_H

#include "SundanceEvaluatableExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      /**
       * LeafExpr is a base class for those evaluatable exprs that are
       * leaves in an expression tree, for example, functions and
       * constants. 
       */
      class LeafExpr : public virtual EvaluatableExpr
        {
        public:
          /** */
          LeafExpr();

          /** */
          virtual ~LeafExpr() {;}

          

          



        protected:
        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
