/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PRODUCTTRANSFORMATION_H
#define SUNDANCE_PRODUCTTRANSFORMATION_H

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
       * ProductTransformation is a base class for any transformation
       * which takes the two operands of a product (left, right) and produces
       * a new expression mathematically equivalent to the original
       * product. This will be used to effect simplification
       * transformations on product expressions.
       */
      class ProductTransformation : public SymbolicTransformation
        {
        public:
          /** */
          ProductTransformation();

          /** */
          virtual ~ProductTransformation(){;}

          /** 
           * Test whether the transform is applicable in this case,
           * and if it is, apply it. The return value is true is the
           * transformation was applied, otherwise false. 
           * Returns by non-const reference
           * the transformed expression. 
           * 
           */
          virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                   RefCountPtr<ScalarExpr>& rtn) const = 0 ;
        };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
