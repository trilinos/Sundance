/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasisFamily.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;




int BasisFamily::order() const 
{
  return ptr()->order();
}

int BasisFamily::dim() const 
{
  return ptr()->dim();
}

bool BasisFamily::operator==(const BasisFamily& other) const
{
  return !(*this < other || other < *this);
}

int BasisFamily::size(const Array<BasisFamily>& b)
{
  int rtn = 0;
  for (int i=0; i<b.size(); i++) rtn += b[i].dim();
  return rtn;
}

int BasisFamily::nNodes(const CellType& cellType) const 
{
  return ptr()->nNodes(cellType);
}

BasisFamily BasisFamily::getBasis(const Expr& expr)
{
  TEST_FOR_EXCEPTION(expr.size() > 1, RuntimeError, "non-scalar expression in BasisFamily::getBasis()");

  /* First try to process assuming the input is an unknown function */
  const UnknownFuncElement* u 
    = dynamic_cast<const UnknownFuncElement*>(expr[0].ptr().get());
  if (u != 0)
    {
      const UnknownFunctionStub* stub = u->master();
      const UnknownFunction* unk = dynamic_cast<const UnknownFunction*>(stub);
      TEST_FOR_EXCEPTION(unk==0, InternalError, 
                         "In BasisFamily::getBasis() input expression " 
                         << expr[0]
                         << " derives from UnknownFunctionStub but not from "
                         "UnknownFunction");
      return unk->basis()[u->myIndex()];
    }

  /* Next try to process assuming the input is a test function */
  const TestFuncElement* t 
    = dynamic_cast<const TestFuncElement*>(expr[0].ptr().get());
  if (t != 0)
    {
      const TestFunctionStub* stub = t->master();
      const TestFunction* test = dynamic_cast<const TestFunction*>(stub);
      TEST_FOR_EXCEPTION(test==0, InternalError, 
                         "In BasisFamily::getBasis() input expression " 
                         << expr[0]
                         << " derives from TestFunctionStub but not from "
                         "TestFunction");
      return test->basis()[t->myIndex()];
    }

  /* Next try to process assuming the input is a discrete function */
  const DiscreteFuncElement* d
    = dynamic_cast<const DiscreteFuncElement*>(expr[0].ptr().get());
  if (d != 0)
    {
      const DiscreteFunctionStub* stub = d->master();
      const DiscreteFunction* disc 
        = dynamic_cast<const DiscreteFunction*>(stub);
      TEST_FOR_EXCEPTION(disc==0, InternalError, 
                         "In BasisFamily::getBasis() input expression " 
                         << expr[0]
                         << " derives from DiscreteFunctionStub but not from "
                         "DiscreteFunction");
      return disc->basis()[d->myIndex()];
    }

  
  
}
