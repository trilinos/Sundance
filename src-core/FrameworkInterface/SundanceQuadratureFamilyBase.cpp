/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;
using namespace TSFExtended;

QuadratureFamilyBase::QuadratureFamilyBase(int order)
: order_(order)
{;}

XMLObject QuadratureFamilyBase::toXML() const
{
  XMLObject rtn("QuadratureFamilyBase");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}
