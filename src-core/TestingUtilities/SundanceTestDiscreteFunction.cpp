/* @HEADER@ */
/* @HEADER@ */


#include "SundanceTestDiscreteFunction.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceTesting;
using namespace TSFExtended;
using namespace Teuchos;
using namespace std;

ADReal TestDiscreteFunction::evaluate() const
{
  return field_.evaluate();
}
