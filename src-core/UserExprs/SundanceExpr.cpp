/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceConstantExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceUnaryMinus.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceStdSumTransformations.hpp"
#include "SundanceStdProductTransformations.hpp"
#include "SundanceNonlinearUnaryOp.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceParameter.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;

Expr::Expr(const double& c)
	: TSFExtended::Handle<ExprBase>(new ConstantExpr(c))
{}



XMLObject Expr::toXML() const
{
  TimeMonitor t(outputTimer());

	return ptr()->toXML();
}

string Expr::toString() const 
{
  TimeMonitor t(outputTimer());

	TeuchosOStringStream ss;
	ptr()->toText(ss, false);
	return TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss);
}

string Expr::toLatex() const 
{
  TimeMonitor t(outputTimer());

  TeuchosOStringStream ss;
	ptr()->toLatex(ss, false);
  //	ss << ends;
	return TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss);
}








Expr Expr::operator+(const Expr& other) const 
{
  TimeMonitor t(opTimer());

  /* if both operands are simple scalars, add them */
  if (size() == 1 && other.size()==1)
    {
      return (*this)[0].sum(other[0], 1);
    }

  /* otherwise, create a list of the sums */
  Array<Expr> rtn(size());
  if (size() == other.size())
    {
      for (int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i] + other[i];
        }
    }
  return new ListExpr(rtn);
}

Expr Expr::operator-(const Expr& other) const 
{
  TimeMonitor t(opTimer());

  /* if both operands are simple scalars, subtract them */
  if (size() == 1 && other.size()==1)
    {
      return (*this)[0].sum(other[0], -1);
    }

  /* otherwise, create a list of the sums */
  Array<Expr> rtn(size());
  if (size() == other.size())
    {
      for (int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i] - other[i];
        }
    }
  return new ListExpr(rtn);
}


Expr Expr::sum(const Expr& other, int sign) const 
{
  RefCountPtr<ScalarExpr> rtn;
  RefCountPtr<ScalarExpr> sThis = rcp_dynamic_cast<ScalarExpr>(ptr());
  RefCountPtr<ScalarExpr> sOther = rcp_dynamic_cast<ScalarExpr>(other.ptr());

  TEST_FOR_EXCEPTION(ptr().get()==NULL, InternalError,
                     "Expr::sum() detected null this pointer");

  TEST_FOR_EXCEPTION(sThis.get()==NULL, InternalError,
                     "Expr::sum(): Left operand " << toString() 
                     << " is a non-scalar expression. All list structure "
                     "should have been handled before this point");

  TEST_FOR_EXCEPTION(sOther.get()==NULL, InternalError,
                     "Expr::sum(): Right operand " << other.toString() 
                     << " is a non-scalar expression. All list structure "
                     "should have been handled before this point");

  static StdSumTransformations trans;

  if (trans.doTransform(sThis, sOther, sign, rtn)) 
    {
      if (SymbolicTransformation::classVerbosity() > 0)
        {
          Out::println("Expr::sum() transformed sum\n[" 
                       + toString() + "+"
                       + other.toString() + "]\n to\n ["
                       + rtn->toString() + "]");
        }
      return handle(rtn);
    }

  else return new SumExpr(sThis, sOther, sign);
}


Expr Expr::operator*(const Expr& other) const 
{
  TimeMonitor t(opTimer());

  /* if both operands are simple scalars, multiply them */
  if (size() == 1 && other.size()==1)
    {
      return (*this)[0].multiply(other[0]);
    }

  /* if the left operand is a scalar, multiply it through */
  if (size()==1)
    {
      Array<Expr> rtn(other.size());
      for (int i=0; i<other.size(); i++)
        {
          rtn[i] = (*this)[0] * other[i];
        }
      return new ListExpr(rtn);
    }

  /* if the right operand is a scalar, multiply it through */
  if (other.size()==1)
    {
      Array<Expr> rtn(size());
      for (int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i] * other[0];
        }
      return new ListExpr(rtn);
    }

  /* if both operands are flat vectors, take their dot product */
  if (size() == totalSize() && other.size()==other.totalSize() 
      && size() == other.size())
    {
      Expr rtn = new ZeroExpr();

      for (int i=0; i<size(); i++)
        {
          rtn = rtn + (*this)[i]*other[i];
        }
      return rtn;
    }

  /* if the left operand is a rectangular matrix and the 
   * right operand is a vector */
  int cols = (*this)[0].size();
  bool rectangular = true;
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i].size() != cols) rectangular = false;
    }
  TEST_FOR_EXCEPTION(!rectangular, BadSymbolicsError,
                     "Expr::operator* detected list-list multiplication "
                     "with a non-rectangular left operator " 
                     << toString());
  
  TEST_FOR_EXCEPTION(cols != other.size(), BadSymbolicsError,
                     "Expr::operator* detected mismatched dimensions in "
                     "list-list multiplication. Left operator is "
                     << toString() << ", right operator is "
                     << other.toString());
  
  Array<Expr> rtn(size());
  for (int i=0; i<size(); i++)
    {
      rtn[i] = (*this)[i] * other;
    }

  return new ListExpr(rtn);
}

Expr Expr::divide(const Expr& other) const 
{
  RefCountPtr<ScalarExpr> sOther = rcp_dynamic_cast<ScalarExpr>(other.ptr());
  Expr recip = new NonlinearUnaryOp(sOther, rcp(new StdReciprocal()));
  return (*this)[0].multiply(recip);
}

Expr Expr::multiply(const Expr& other) const 
{
  RefCountPtr<ScalarExpr> rtn;
  RefCountPtr<ScalarExpr> sThis = rcp_dynamic_cast<ScalarExpr>(ptr());
  RefCountPtr<ScalarExpr> sOther = rcp_dynamic_cast<ScalarExpr>(other.ptr());

  TEST_FOR_EXCEPTION(sThis.get()==NULL, InternalError,
                     "Expr::multiply(): Left operand " << toString() 
                     << " is a non-scalar expression. All list structure "
                     "should have been handled before this point");

  TEST_FOR_EXCEPTION(sOther.get()==NULL, InternalError,
                     "Expr::multiply(): Right operand " << other.toString() 
                     << " is a non-scalar expression. All list structure "
                     "should have been handled before this point");

  static StdProductTransformations trans;

  if (trans.doTransform(sThis, sOther, rtn)) 
    {
      if (SymbolicTransformation::classVerbosity() > 0)
        {
          Out::println("Expr::operator*() transformed product\n[" 
                       + toString() + "*"
                       + other.toString() + "]\n to\n ["
                       + rtn->toString() + "]");
        }
      return handle(rtn);
    }

  return new ProductExpr(sThis, sOther);
}

Expr Expr::operator-() const 
{
  TimeMonitor t(opTimer());

  /* if we are a scalar, form a unary minus */
  if (size()==1)
    {
      RefCountPtr<ScalarExpr> sThis 
        = rcp_dynamic_cast<ScalarExpr>((*this)[0].ptr());
      TEST_FOR_EXCEPTION(sThis.get()==NULL, InternalError,
                         "Expr::operator-(): Operand " << (*this)[0].toString() 
                         << " is a non-scalar expression. All list structure "
                         "should have been handled before this point");
      return new UnaryMinus(sThis);
    }

  /* otherwise, distribute the sign change over the list */
  Array<Expr> rtn(size());
  for (int i=0; i<size(); i++)
    {
      rtn[i] = -((*this)[i]);
    }
  return new ListExpr(rtn);
}


Expr Expr::operator/(const Expr& other) const 
{
  TimeMonitor t(opTimer());

  /* if the right operand is a list, this operation
   * makes no sense */
  TEST_FOR_EXCEPTION(other.size() != 1,
                     BadSymbolicsError, 
                     "Expr::operator/ detected division by a non-scalar "
                     "expression " << toString());

  /* if we are a scalar, do simple scalar division */
  if (size()==1)
    {
      return (*this)[0].divide(other[0]);
    }

  /* otherwise, divide each element of the left by the right operand */
  Array<Expr> rtn(size());
  for (int i=0; i<size(); i++)
    {
      rtn[i] = (*this)[i] / other;
    }
  return new ListExpr(rtn);
}

const Expr& Expr::operator[](int i) const
{
  TEST_FOR_EXCEPTION(ptr().get()==NULL, InternalError,
                     "null this detected in Expr::operator[].");

  TEST_FOR_EXCEPTION(i<0 || i>=size(), RuntimeError,
                     "Expr::operator[]() index i=" << i << " out of range [0, "
                     << size() << " in expr " << toString());

  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      const Expr& rtn = le->element(i);
      if (rtn.size()==1) 
        {
          TEST_FOR_EXCEPTION(rtn[0].ptr().get()==NULL, InternalError,
                             "null return detected in Expr::operator[]. This="
                             << toString() << ", i=" << i);
          return rtn[0];
        }
      TEST_FOR_EXCEPTION(rtn.ptr().get()==NULL, InternalError,
                         "null return detected in Expr::operator[]. This="
                         << toString() << ", i=" << i);
      return rtn;
    }
  return *this;
}

void Expr::append(const Expr& expr)
{
  ListExpr* le = dynamic_cast<ListExpr*>(ptr().get());

  if (le != 0)
    {
      le->append(expr);
      return;
    }
  else
    {
      if (ptr().get()==0)
        {
          Array<Expr> e(1);
          e[0] = expr;
          ptr() = rcp(new ListExpr(e));
        }
      else
        {
          Array<Expr> e(2);
          e[0] = *this;
          e[1] = expr;
          ptr() = rcp(new ListExpr(e));
        }
    }
}

Expr Expr::flatten() const
{
  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      return le->flatten();
    }
  return *this;
}

int Expr::size() const
{
  if (ptr().get()==0) return 0;

  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      return le->size();
    }
  return 1;
}

int Expr::totalSize() const
{
  if (ptr().get()==0) return 0;

  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      return le->totalSize();
    }
  return 1;
}

void Expr::setParameterValue(const double& value)
{
  Parameter* pe = dynamic_cast<Parameter*>(ptr().get());
  TEST_FOR_EXCEPTION(pe==0, RuntimeError, 
                     "Expr " << *this << " is not a Parameter expr, and "
                     "so setParameterValue() should not be called");
  pe->setValue(value);
}

namespace SundanceCore
{
  using namespace SundanceUtils;
  Expr List(const Expr& a)
  {
    return new ListExpr(tuple(a));
  }

  Expr List(const Expr& a, const Expr& b)
  {
    return new ListExpr(tuple(a,b));
  }

  Expr List(const Expr& a, const Expr& b, const Expr& c)
  {
    return new ListExpr(tuple(a,b,c));
  }

  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d)
  {
    return new ListExpr(tuple(a,b,c,d));
  }

  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d, const Expr& e)
  {
    return new ListExpr(tuple(a,b,c,d,e));
  }

  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d, const Expr& e, const Expr& f)
  {
    return new ListExpr(tuple(a,b,c,d,e,f));
  }
}

