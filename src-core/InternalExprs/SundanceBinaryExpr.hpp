/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BINARYEXPR_H
#define SUNDANCE_BINARYEXPR_H

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
      /** 
       * BinaryExpr is a base class for binary expressions, e.g., sums
       * and products. It provides a number of helper methods.
       */
      class BinaryExpr : public ExprWithChildren
        {
        public:
          /** construct with left and right operands */
          BinaryExpr(const RefCountPtr<ScalarExpr>& left,
                     const RefCountPtr<ScalarExpr>& right, int sign);

          /** virtual dtor */
          virtual ~BinaryExpr() {;}

          /** */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;

          /** */
          Expr left() const {return child(0);}

          /** */
          Expr right() const {return child(1);}

          /** */
          int sign() const {return sign_;}

          /** Downcast the left expr to an evaluatable expr */
          const EvaluatableExpr* leftEvaluatable() const 
          {return evaluatableChild(0);}

          /** Downcast the right expr to an evaluatable expr */
          const EvaluatableExpr* rightEvaluatable() const 
          {return evaluatableChild(1);}

          /** Downcast the left expr to a scalar expr */
          const ScalarExpr* leftScalar() const {return scalarChild(0);}

          /** Downcast the right expr to a scalar expr */
          const ScalarExpr* rightScalar() const {return scalarChild(1);}

        protected:

          

          /** */
          virtual bool parenthesizeSelf() const = 0 ;
          /** */
          virtual bool parenthesizeOperands() const = 0 ;
          /** */
          virtual const string& xmlTag() const = 0 ;
          /** */
          virtual const string& opChar() const = 0 ;



        private:
          

          int sign_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
