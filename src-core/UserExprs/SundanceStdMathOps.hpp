/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDMATHOPS_H
#define SUNDANCE_STDMATHOPS_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceStdMathFunctors.hpp"
#include "SundanceNonlinearUnaryOp.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY


using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

#define SUNDANCE_UNARY_OP(opName, functorName, description) \
/** \relates Expr description */\
inline Expr opName(const Expr& expr) \
{\
RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());\
    TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,\
                       "non-scalar argument in " #opName " function");\
    return new NonlinearUnaryOp(arg, rcp(new functorName()));\
}

namespace SundanceCore
{
  inline Expr pow(const Expr& expr, const double& p)
  {
    RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
    TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,
                       "non-scalar argument in pow function");
    return new NonlinearUnaryOp(arg, rcp(new PowerFunctor(p)));
  }

  /** \name Elementary math functions */
  //@{
  SUNDANCE_UNARY_OP(reciprocal, StdReciprocal, "reciprocal function");

  SUNDANCE_UNARY_OP(fabs, StdFabs, "absolute value");

  SUNDANCE_UNARY_OP(sign, StdSign, "sign function");

  SUNDANCE_UNARY_OP(exp, StdExp, "exponential function");

  SUNDANCE_UNARY_OP(log, StdLog, "logarithm");

  SUNDANCE_UNARY_OP(sqrt, StdSqrt, "square root"]);

  SUNDANCE_UNARY_OP(sin, StdSin, "sine function");

  SUNDANCE_UNARY_OP(cos, StdCos, "cosine function");

  SUNDANCE_UNARY_OP(tan, StdTan, "tangent function");

  SUNDANCE_UNARY_OP(asin, StdASin, "inverse sine");

  SUNDANCE_UNARY_OP(acos, StdACos, "inverse cosine");

  SUNDANCE_UNARY_OP(atan, StdATan, "inverse tangent");

  SUNDANCE_UNARY_OP(sinh, StdSinh, "hyperbolic sine");

  SUNDANCE_UNARY_OP(cosh, StdCosh, "hyperbolic cosine");

  SUNDANCE_UNARY_OP(tanh, StdTanh, "hyperbolic tangent");

  SUNDANCE_UNARY_OP(asinh, StdASinh, "inverse hyperbolic sine");

  SUNDANCE_UNARY_OP(acosh, StdACosh, "inverse hyperbolic cosine");

  SUNDANCE_UNARY_OP(atanh, StdATanh, "inverse hyperbolic tangent");
//@}

}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
