/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;

SymbolicTransformation::SymbolicTransformation()
{}

RefCountPtr<ScalarExpr> SymbolicTransformation::chooseSign(int sign, 
                                                           const RefCountPtr<ScalarExpr>& expr) 
{
  /* return expr if sign == 1, -expr if sign == -1. No other
   * cases should happen. */
  switch(sign)
    {
    case 1:
      return expr;
    case -1:
      {
        Expr e = -Expr::handle(expr);
        RefCountPtr<ScalarExpr> rtn = rcp_dynamic_cast<ScalarExpr>(e.ptr());
        TEST_FOR_EXCEPTION(rtn.get() == NULL, InternalError,
                           "Non-scalar expr "
                           << e.toString() 
                           << " detected in SymbolicTransformation::chooseSign");
        return rtn;
      }
    default:
      TEST_FOR_EXCEPTION(true, InternalError, 
                         "sign != +/- 1 in Expr::transformSign()");
    }
  return expr;
}

Expr SymbolicTransformation::chooseSign(int sign, 
                                        const Expr& expr) 
{
  /* return expr if sign == 1, -expr if sign == -1. No other
   * cases should happen. */
  switch(sign)
    {
    case 1:
      return expr;
    case -1:
      return -expr;
    default:
      TEST_FOR_EXCEPTION(true, InternalError, 
                         "sign != +/- 1 in Expr::transformSign()");
    }
  return expr;
}

RefCountPtr<ScalarExpr> SymbolicTransformation::getScalar(const Expr& expr)
{
  RefCountPtr<ScalarExpr> s = rcp_dynamic_cast<ScalarExpr>(expr.ptr());

  TEST_FOR_EXCEPTION(s.get()==NULL, InternalError,
                     "non-scalar detected in SymbolicTransformation::getScalar");

  return s;
}
