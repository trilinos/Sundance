/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNARYEXPR_H
#define SUNDANCE_UNARYEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprWithChildren.hpp"
#include "SundanceExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** */
      class UnaryExpr : public ExprWithChildren
        {
        public:
          /** construct with the argument */
          UnaryExpr(const RefCountPtr<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryExpr() {;}

          /** Return the operand */
          Expr arg() const {return child(0);}

          /** Downcast the argument to an evaluatable expr */
          const EvaluatableExpr* evaluatableArg() const 
          {return evaluatableChild(0);}

          /** Indicate whether the given derivative of this expression
           * is nonzero. For all non-DiffOp unary operators, the
           * derivative is nonzero if the operand's derivative is 
           * nonzero. DiffOp will have a specialized implementation of 
           * this method. */
          virtual bool hasNonzeroDeriv(const MultipleDeriv& f) const
            {return evaluatableArg()->hasNonzeroDeriv(f);}

          

          /** Return the index by which the given deriv set is known
           * by the operand. */
          int argDerivSetIndex(int derivSetIndex) const 
          {return childDerivSetIndex(0, derivSetIndex);}

        protected:
        private:
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
