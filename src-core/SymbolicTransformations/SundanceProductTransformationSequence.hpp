/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PRODUCTTRANSFORMATIONSEQUENCE_H
#define SUNDANCE_PRODUCTTRANSFORMATIONSEQUENCE_H

#include "SundanceDefs.hpp"
#include "SundanceProductTransformation.hpp"
#include "Teuchos_Array.hpp"

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
       * ProductTransformationSequence is a sequence of transformations
       * to be applied to a product, producing a transformed expression. 
       */
      class ProductTransformationSequence : public ProductTransformation,
                                        public Array<RefCountPtr<ProductTransformation> >
        {
        public:
          /** */
          ProductTransformationSequence();

          /** */
          virtual ~ProductTransformationSequence(){;}

          /**
           * Test whether the transform is applicable in this case,
           * and if it is, apply it. The return value is true is the
           * transformation was applied, otherwise false. 
           * Returns by non-const reference
           * the transformed expression. 
           *
           * For ProductTransformationSequence, this is implemented by
           * trying to apply all transformations in sequence. If one
           * succeeds, we exit immediately with true.
           */
          virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                   RefCountPtr<ScalarExpr>& rtn) const ;

        
        };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
