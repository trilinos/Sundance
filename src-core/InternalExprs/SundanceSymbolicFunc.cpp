/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSymbolicFunc.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

SymbolicFunc::SymbolicFunc()
  : ListExpr()
{}


void SymbolicFunc::substituteZero() const 
{
  for (int i=0; i<size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEST_FOR_EXCEPTION(u==0, InternalError, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteZero()");
      u->substituteZero();
    }
}

void SymbolicFunc
::substituteFunction(const RefCountPtr<DiscreteFunctionStub>& u0) const
{
  TEST_FOR_EXCEPTION(size() != u0->size(), InternalError,
                     "Mismatch between sizes of symbolic " << toString()
                     << " and discrete func " << u0->toString()
                     << " in substituteFunction()");

  for (int i=0; i<size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEST_FOR_EXCEPTION(u==0, InternalError, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteFunction()");

      RefCountPtr<DiscreteFuncElement> df 
        = rcp_dynamic_cast<DiscreteFuncElement>(u0->element(i).ptr());
      u->substituteFunction(df);
    }
}

