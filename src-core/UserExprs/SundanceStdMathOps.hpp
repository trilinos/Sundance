/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDMATHOPS_H
#define SUNDANCE_STDMATHOPS_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceNonlinearUnaryOp.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY


using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

#define DEFINE_UNARY_OP(opName, functorName, description, funcDefinition, derivDefinition) \
  /** */\
class functorName : public Internal::UnaryFunctor \
{\
public:\
/** ctor for description functor */\
functorName() : Internal::UnaryFunctor(#opName) {;} \
/** Evaluate function at an array of values */ \
void eval(const double* const x, int nx, double* f) const ;\
/** Evaluate function and derivative at an array of values */\
void eval(const double* const x, \
                  int nx, \
                  double* f, \
                  double* df) const ;\
};\
/** \relates Expr description */\
inline Expr opName(const Expr& expr) \
{\
RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());\
    TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,\
                       "non-scalar argument in " #opName " function");\
    return new NonlinearUnaryOp(arg, rcp(new functorName()));\
}\
inline void functorName::eval(const double* const x, int nx, double* f) const \
{\
  for (int i=0; i<nx; i++) f[i] = funcDefinition;\
}\
inline void functorName::eval(const double* const x, \
                  int nx, \
                  double* f, \
                  double* df) const \
{ \
  for (int i=0; i<nx; i++) \
    { \
      f[i] = funcDefinition;\
      df[i] = derivDefinition;\
    }\
} 



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  /** */
  class PowerFunctor : public Internal::UnaryFunctor
  {
  public:
    /** */
    PowerFunctor(const double& p);
    
    /** Evaluate power function and deriv at an array of values */ 
    void eval(const double* const x, 
                  int nx, 
                  double* f, 
                  double* df) const ;
    /** Evaluate power function at an array of values */ 
    void eval(const double* const x, int nx, double* f) const ;
  private:
    double p_;
  };

  inline Expr pow(const Expr& expr, const double& p)
  {
    RefCountPtr<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
    TEST_FOR_EXCEPTION(arg.get()==0, RuntimeError,
                       "non-scalar argument in pow function");
    return new NonlinearUnaryOp(arg, rcp(new PowerFunctor(p)));
  }

  using std::string;
  using std::ostream;

  DEFINE_UNARY_OP(reciprocal, StdReciprocal, "reciprocal function", 1.0/x[i], -f[i]*f[i]);

  DEFINE_UNARY_OP(exp, StdExp, "exponential function", ::exp(x[i]), f[i]);

  DEFINE_UNARY_OP(log, StdLog, "logarithm", ::log(x[i]), 1.0/x[i]);

  DEFINE_UNARY_OP(sqrt, StdSqrt, "square root", ::sqrt(x[i]), 0.5/f[i]);

  DEFINE_UNARY_OP(sin, StdSin, "sine function", ::sin(x[i]), ::cos(x[i]));

  DEFINE_UNARY_OP(cos, StdCos, "cosine function", ::cos(x[i]), -::sin(x[i]));

  DEFINE_UNARY_OP(tan, StdTan, "tangent function", 
                  ::tan(x[i]), 1.0 + f[i]*f[i]);

  DEFINE_UNARY_OP(asin, StdASin, "inverse sine", 
                  ::asin(x[i]), 1.0/::sqrt(1.0-x[i]*x[i]));

  DEFINE_UNARY_OP(acos, StdACos, "inverse cosine", 
                  ::acos(x[i]), -1.0/::sqrt(1.0-x[i]*x[i]));

  DEFINE_UNARY_OP(atan, StdATan, "inverse tangent", 
                  ::atan(x[i]), 1.0/(1.0 + x[i]*x[i]));

  DEFINE_UNARY_OP(sinh, StdSinh, "hyperbolic sine",
                  ::sinh(x[i]), ::cosh(x[i]));

  DEFINE_UNARY_OP(cosh, StdCosh, "hyperbolic cosine",
                  ::cosh(x[i]), ::sinh(x[i]));

  DEFINE_UNARY_OP(tanh, StdTanh, "hyperbolic tangent",
                  ::tanh(x[i]), 1.0 - f[i]*f[i]);

  DEFINE_UNARY_OP(asinh, StdASinh, "inverse hyperbolic sine",
                  ::asinh(x[i]), 1.0/::sqrt(1.0 + x[i]*x[i]));

  DEFINE_UNARY_OP(acosh, StdACosh, "inverse hyperbolic cosine",
                  ::acosh(x[i]), 1.0/::sqrt(x[i]*x[i]-1.0));

  DEFINE_UNARY_OP(atanh, StdATanh, "inverse hyperbolic tangent",
                  ::atanh(x[i]), 1.0/(1.0 - x[i]*x[i]));


}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
