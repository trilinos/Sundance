/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFunctionBase.hpp"
#include "SundanceDiscreteFuncElement.hpp"


using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

DiscreteFunctionBase::DiscreteFunctionBase(const string& name, int nElems)
	: ListExpr()
{
  for (int i=0; i<nElems; i++)
    {
      string elemName = name + "[" + Teuchos::toString(i) + "]";
      append(new DiscreteFuncElement(this, elemName, i));
    }
}


