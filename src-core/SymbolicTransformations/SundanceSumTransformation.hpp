/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMTRANSFORMATION_H
#define SUNDANCE_SUMTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceSymbolicTransformation.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY 

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** 
       * SumTransformation is a base class for any transformation
       * which takes the two operands of a sum (left, right) and produces
       * a new expression mathematically equivalent to the original
       * sum. This will be used to effect simplification
       * transformations on sum expressions.
       */
      class SumTransformation : public SymbolicTransformation
        {
        public:
          /** */
          SumTransformation();

          /** */
          virtual ~SumTransformation(){;}

          /** 
           * Test whether the transform is applicable in this case,
           * and if it is, apply it. The return value is true is the
           * transformation was applied, otherwise false. 
           * Returns by non-const reference
           * the transformed expression. 
           * 
           */
          virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                   int sign, RefCountPtr<ScalarExpr>& rtn) const = 0 ;
        };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
