/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceSumTransformationSequence.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;

SumTransformationSequence::SumTransformationSequence()
  : SumTransformation(), 
    Array<RefCountPtr<SumTransformation> >()
{;}

bool SumTransformationSequence::doTransform(const RefCountPtr<ScalarExpr>& left, 
                                            const RefCountPtr<ScalarExpr>& right,
                                            int sign, RefCountPtr<ScalarExpr>& rtn) const
{
  for (unsigned int i=0; i<size(); i++)
    {
      if ((*this)[i]->doTransform(left, right, sign, rtn)) return true;
    }

  return false;
}

