/* @HEADER@ */
/* @HEADER@ */

#include "SundanceGaussianQuadrature.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

GaussianQuadrature::GaussianQuadrature(int order)
  : QuadratureFamilyStub(order)
{;}

XMLObject GaussianQuadrature::toXML() const 
{
  XMLObject rtn("GaussianQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}
