/* @HEADER@ */
/* @HEADER@ */

#include "SundanceLeafExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

LeafExpr::LeafExpr()
  : EvaluatableExpr()
{}

int LeafExpr::setupEval(const RegionQuadCombo& region,
                        const EvaluatorFactory* factory) const
{
  /* first check to see if rule tables for this set of derivatives 
   * have already been created. If so, we're done here. */
  if (checkForKnownRegion(region))
    {
      return getDerivSetIndex(region);
    }
  
  /* OK, if we've not returned, we need to create evaluation rules 
   * for this deriv set. */
  int derivSetIndex = registerRegion(region, currentDerivSuperset(), 
                                      factory);

  return derivSetIndex;
}




