/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBoundaryCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool BoundaryCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const BoundaryCellPredicate*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to BoundaryCellPredicate::lessThan() should be "
                     "a BoundaryCellPredicate pointer.");

  return false;
}

bool BoundaryCellPredicate::test(int cellLID) const 
{
  return mesh().numCofacets(cellDim(), cellLID) == 1;
}

XMLObject BoundaryCellPredicate::toXML() const 
{
  XMLObject rtn("BoundaryCellPredicate");
  return rtn;
}

