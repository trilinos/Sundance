/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNullCellFilterBase.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;
using namespace TSFExtended;

NullCellFilterBase::NullCellFilterBase()
{;}

bool NullCellFilterBase::lessThan(const CellFilterBase* other) const
{
  const NullCellFilterBase* ncf = dynamic_cast<const NullCellFilterBase*>(other);
  TEST_FOR_EXCEPTION(ncf==0, RuntimeError,
                     "argument " << other->describe() 
                     << " to NullCellFilter::lessThan() could not be cast "
                     "to a NullCellFilter*");
  /* All null cell filters are equivalent */
  return false;
}


XMLObject NullCellFilterBase::toXML() const
{
  XMLObject rtn("NullCellFilterBase");
  return rtn;
}
