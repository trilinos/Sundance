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
       * constants. It exists in order to provide the non-recursive
       * implementation of setupEval() appropriate to leaves.
       */
      class LeafExpr : public virtual EvaluatableExpr
        {
        public:
          /** */
          LeafExpr();

          /** */
          virtual ~LeafExpr() {;}

          /** Create an evaluator for this expression. Since this is
           * a leaf, this method is non-recursive, i.e., it does
           * not have any operands to call setupEval() on. */
          virtual int setupEval(const EvalContext& region,
                                const EvaluatorFactory* factory,
                                bool regardFuncAsConstant) const ;


        protected:
        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
