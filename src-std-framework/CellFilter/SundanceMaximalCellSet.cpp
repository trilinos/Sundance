/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMaximalCellSet.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceImplicitCellSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

MaximalCellSet::MaximalCellSet()
  : CellSetBase()
{;}

XMLObject MaximalCellSet::toXML() const 
{
  XMLObject rtn("MaximalCellSet");
  return rtn;
}

bool MaximalCellSet::lessThan(const CellSetBase* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const MaximalCellSet*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to MaximalCellSet::lessThan() should be "
                     "a MaximalCellSet pointer.");

  return false;
}
