/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellFilterBase.hpp"
#include "SundanceNullCellFilterBase.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;
using namespace TSFExtended;

CellFilterBase::CellFilterBase()
{;}

XMLObject CellFilterBase::toXML() const
{
  XMLObject rtn("CellFilterBase");
  rtn.addAttribute("id", Teuchos::toString(id()));
  return rtn;
}

RefCountPtr<CellFilterBase> CellFilterBase::makeNullDomain() const
{
  return rcp(new NullCellFilterBase());
}

bool CellFilterBase::isNullDomain() const
{
  return 0 != dynamic_cast<const NullCellFilterBase*>(this);
}
