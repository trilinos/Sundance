/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceProductTransformationSequence.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;

ProductTransformationSequence::ProductTransformationSequence()
  : ProductTransformation(), 
    Array<RefCountPtr<ProductTransformation> >()
{;}

bool ProductTransformationSequence::doTransform(const RefCountPtr<ScalarExpr>& left, 
                                                const RefCountPtr<ScalarExpr>& right,
                                                RefCountPtr<ScalarExpr>& rtn) const
{
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]->doTransform(left, right, rtn)) return true;
    }

  return false;
}

