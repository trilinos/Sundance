/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureFamily.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



XMLObject QuadratureFamily::toXML() const 
{
  return ptr()->toXML();
}

int QuadratureFamily::order() const 
{
  return ptr()->order();
}
