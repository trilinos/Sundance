/* @HEADER@ */
/* @HEADER@ */


#include "SundanceTestUnknownFunction.hpp"
#include "SundanceTestDiscreteFunction.hpp"
#include "SundanceUnknownFuncElement.hpp"
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

Expr TestUnknownFunction::createDiscreteFunction() const 
{
  const UnknownFuncElement* u 
    = dynamic_cast<const UnknownFuncElement*>(element(0).ptr().get());

  return new TestDiscreteFunction(field_, u->rootName() + "0");
}
