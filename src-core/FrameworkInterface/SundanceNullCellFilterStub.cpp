/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNullCellFilterStub.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

NullCellFilterStub::NullCellFilterStub()
{;}

bool NullCellFilterStub::lessThan(const CellFilterStub* other) const
{
  const NullCellFilterStub* ncf = dynamic_cast<const NullCellFilterStub*>(other);
  TEST_FOR_EXCEPTION(ncf==0, RuntimeError,
                     "argument " << other->describe() 
                     << " to NullCellFilter::lessThan() could not be cast "
                     "to a NullCellFilter*");
  /* All null cell filters are equivalent */
  return false;
}


XMLObject NullCellFilterStub::toXML() const
{
  XMLObject rtn("NullCellFilterStub");
  return rtn;
}
