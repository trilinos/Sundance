/* @HEADER@ */
/* @HEADER@ */

#include "SundanceStdMathOps.hpp"
#include "SundanceNonlinearUnaryOp.hpp"
#include "SundanceScalarExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


namespace SundanceCore
{
Expr sin(const Expr& expr)
{
  RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
  TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,
                     "non-scalar argument in sin() function");
  return new NonlinearUnaryOp(arg, rcp(new StdSin()));
}

Expr cos(const Expr& expr)
{
  RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
  TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,
                     "non-scalar argument in cos() function");
  return new NonlinearUnaryOp(arg, rcp(new StdCos()));
}

Expr exp(const Expr& expr)
{
  RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
  TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,
                     "non-scalar argument in exp() function");
  return new NonlinearUnaryOp(arg, rcp(new StdExp()));
}
}






