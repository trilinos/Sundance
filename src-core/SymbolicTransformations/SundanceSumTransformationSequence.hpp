/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUMTRANSFORMATIONSEQUENCE_H
#define SUNDANCE_SUMTRANSFORMATIONSEQUENCE_H

#include "SundanceDefs.hpp"
#include "SundanceSumTransformation.hpp"
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
       * SumTransformationSequence is a sequence of transformations
       * to be applied to a sum, producing a transformed expression. 
       */
      class SumTransformationSequence : public SumTransformation,
                                        public Array<RefCountPtr<SumTransformation> >
        {
        public:
          /** */
          SumTransformationSequence();

          /** */
          virtual ~SumTransformationSequence(){;}

          /**
           * Test whether the transform is applicable in this case,
           * and if it is, apply it. The return value is true is the
           * transformation was applied, otherwise false. 
           * Returns by non-const reference
           * the transformed expression. 
           *
           * For SumTransformationSequence, this is implemented by
           * trying to apply all transformations in sequence. If one
           * succeeds, we exit immediately with true.
           */
          virtual bool doTransform(const RefCountPtr<ScalarExpr>& left, const RefCountPtr<ScalarExpr>& right,
                                   int sign, RefCountPtr<ScalarExpr>& rtn) const ;

        
        };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
