/* @HEADER@ */
/* @HEADER@ */

#include "SundanceTestFunctionBase.hpp"
#include "SundanceTestFuncElement.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

TestFunctionBase::TestFunctionBase(const string& name, int nElems)
	: SymbolicFunc()
{
  for (int i=0; i<nElems; i++)
    {
      string elemName = name + "[" + Teuchos::toString(i) + "]";
      append(new TestFuncElement(this, elemName, i));
    }
}



