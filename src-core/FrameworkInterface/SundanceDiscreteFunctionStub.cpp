/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"


using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteFunctionStub::DiscreteFunctionStub(const string& name, int nElems)
	: ListExpr()
{
  for (int i=0; i<nElems; i++)
    {
      string suffix = "[" + Teuchos::toString(i) + "]";
      append(new DiscreteFuncElement(this, name, suffix, i));
    }
}


