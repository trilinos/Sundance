/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMEXPR_H
#define SUNDANCE_SUMEXPR_H

#include "SundanceBinaryExpr.hpp"
#include "SundanceSumEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /**
       * SumExpr is the internal representation of an addition or subtraction
       * node in the expression tree.
       */
      class SumExpr : public BinaryExpr,
                      public GenericEvaluatorFactory<SumExpr, SumEvaluator>
        {
        public:
          /** */
          SumExpr(const RefCountPtr<ScalarExpr>& a, 
                  const RefCountPtr<ScalarExpr>& b, int sign);

          /** virtual dtor */
          virtual ~SumExpr() {;}

          /** In order to return true, both left and right operands
           * must contain test functions.  */
          virtual bool allTermsHaveTestFunctions() const ;

          /** */
          virtual bool isHungryDiffOp() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        protected:
          /** */
          virtual bool parenthesizeSelf() const {return true;}
          /** */
          virtual bool parenthesizeOperands() const {return false;}
          /** */
          virtual const string& xmlTag() const ;
          /** */
          virtual const string& opChar() const ;

        private:


        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
