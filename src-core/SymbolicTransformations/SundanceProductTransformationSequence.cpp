/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceProductTransformationSequence.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceOut.hpp"

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
  SUNDANCE_OUT(verbosity() > VerbMedium,
               "testing whether to transform product: " << endl
               << "left = " << left->toString() << endl
               << "right = " << right->toString());
  
  for (unsigned int i=0; i<size(); i++)
    {
      if ((*this)[i]->doTransform(left, right, rtn)) return true;
    }

  return false;
}

