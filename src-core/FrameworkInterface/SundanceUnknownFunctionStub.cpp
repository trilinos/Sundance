/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;



UnknownFunctionStub::UnknownFunctionStub(const string& name, int nElems)
	: SymbolicFunc()
{
  for (int i=0; i<nElems; i++)
    {
      string suffix = "[" + Teuchos::toString(i) + "]";
      append(new UnknownFuncElement(this, name, suffix, i));
    }
}



