/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnknownFunctionBase.hpp"
#include "SundanceUnknownFuncElement.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;



UnknownFunctionBase::UnknownFunctionBase(const string& name, int nElems)
	: SymbolicFunc()
{
  for (int i=0; i<nElems; i++)
    {
      string myName = name + "[" + Teuchos::toString(i) + "]";
      append(new UnknownFuncElement(this,  myName, i));
    }
}



