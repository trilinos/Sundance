/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilterStub.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellFilterStub::CellFilterStub()
{;}

XMLObject CellFilterStub::toXML() const
{
  XMLObject rtn("CellFilterStub");
  rtn.addAttribute("id", Teuchos::toString(id()));
  return rtn;
}

RefCountPtr<CellFilterStub> CellFilterStub::makeNullRegion() const
{
  return rcp(new NullCellFilterStub());
}

bool CellFilterStub::isNullRegion() const
{
  return 0 != dynamic_cast<const NullCellFilterStub*>(this);
}
