/* @HEADER@ */
/* @HEADER@ */

#include "SundancePositionalCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool PositionalCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const PositionalCellPredicate*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to PositionalCellPredicate::lessThan() should be "
                     "a PositionalCellPredicate pointer.");

  return func_ < dynamic_cast<const PositionalCellPredicate*>(other)->func_;
}

bool PositionalCellPredicate::test(int cellLID) const 
{
  if (cellDim()==0)
    {
      return (*func_)(mesh().nodePosition(cellLID));
    }
  else
    {
      Array<int> facets;
      mesh().getFacetArray(cellDim(), cellLID, 0, facets);

      for (int i=0; i<facets.size(); i++)
        {
          if ((*func_)(mesh().nodePosition(facets[i])) == false) return false;
        }
      return true;
    }
}

XMLObject PositionalCellPredicate::toXML() const 
{
  XMLObject rtn("PositionalCellPredicate");
  return rtn;
}

