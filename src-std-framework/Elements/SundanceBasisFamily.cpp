/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasisFamily.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



XMLObject BasisFamily::toXML() const 
{
  return ptr()->toXML();
}

int BasisFamily::order() const 
{
  return ptr()->order();
}
