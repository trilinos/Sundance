/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceOut.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

QuadratureFamilyStub::QuadratureFamilyStub(int order)
: order_(order)
{;}

XMLObject QuadratureFamilyStub::toXML() const
{
  XMLObject rtn("QuadratureFamilyStub");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}
