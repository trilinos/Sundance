/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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
#include "SundanceComplexExpr.hpp"
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

Expr::Expr(const complex<double>& c)
	: TSFExtended::Handle<ExprBase>(new ComplexExpr(new ConstantExpr(c.real()),
                                                  new ConstantExpr(c.imag())))
{}

bool Expr::isComplex() const
{
  return dynamic_cast<const ComplexExpr*>(ptr().get()) != 0;
}

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

bool Expr::sameAs(const Expr& other) const 
{
  if (this->lessThan(other)) return false;
  if (other.lessThan(*this)) return false;
  return true;
}

bool Expr::operator<(const Expr& other) const 
{
  return this->lessThan(other);
}

bool Expr::lessThan(const Expr& other) const
{
  RefCountPtr<ScalarExpr> sThis = rcp_dynamic_cast<ScalarExpr>(ptr());
  RefCountPtr<ScalarExpr> sOther = rcp_dynamic_cast<ScalarExpr>(other.ptr());

  TEST_FOR_EXCEPTION(sThis.get()==0, InternalError,
                     "ordering not defined for non-scalar expression "
                     << toString());

  TEST_FOR_EXCEPTION(sOther.get()==0, InternalError,
                     "ordering not defined for non-scalar expressions"
                     << other.toString());

  
  const ConstantExpr* cMe = dynamic_cast<const ConstantExpr*>(sThis.get());
  const ConstantExpr* cOther = dynamic_cast<const ConstantExpr*>(sOther.get());

  /* constants are placed before any other expr type */
  if (cMe != 0 && cOther==0) return true;
  if (cOther != 0 && cMe==0) return false;
  if (cOther != 0 && cMe != 0) return cMe->lessThan(cOther);

  /* Move generic spatial constants, e.g., parameters to the left.
   * Because the values might change with time, we can't order on values.  */
  const SpatiallyConstantExpr* scMe 
    = dynamic_cast<const SpatiallyConstantExpr*>(sThis.get());
  const SpatiallyConstantExpr* scOther 
    = dynamic_cast<const SpatiallyConstantExpr*>(sOther.get());
  if (scMe != 0 && scOther==0) return true;
  if (scOther != 0 && scMe==0) return false;


  /* try ordering non-constant exprs by type name */
  if (sThis->typeName() < sOther->typeName()) return true;
  if (sThis->typeName() > sOther->typeName()) return false;

  /* if type names are the same, go to base class to do comparison */
  return sThis->lessThan(sOther.get());
}

SundanceUtils::Map<Expr, int> Expr::getSumTree() const
{
  RefCountPtr<ScalarExpr> sThis = rcp_dynamic_cast<ScalarExpr>(ptr());

  TEST_FOR_EXCEPTION(sThis.get()==0, InternalError,
                     "etSumTree() not defined for non-scalar expression "
                     << toString());
  
  const SumExpr* s = dynamic_cast<const SumExpr*>(sThis.get());
  const UnaryMinus* u = dynamic_cast<const UnaryMinus*>(sThis.get());
  if (s != 0)
    {
      return s->getSumTree();
    }
  else if (u != 0)
    {
      SundanceUtils::Map<Expr, int> rtn;
      rtn.put(u->arg(), -1);
      return rtn;
    }
  else
    {
      SundanceUtils::Map<Expr, int> rtn;
      rtn.put(*this, 1);
      return rtn;
    }
  
}

Expr Expr::operator+(const Expr& other) const 
{
  TimeMonitor t(opTimer());

  /* if both operands are simple scalars, add them */
  if (this->size() == 1 && other.size()==1)
    {
      if (!(*this)[0].isComplex() && !other[0].isComplex())
        {
          return (*this)[0].sum(other[0], 1);
        }
      else 
        {
          Expr rtn = new ComplexExpr((*this)[0].real() + other[0].real(),
                                     (*this)[0].imag() + other[0].imag());
          return rtn;
        }
    }

  /* otherwise, create a list of the sums */
  Array<Expr> rtn(this->size());
  if (this->size() == other.size())
    {
      for (unsigned int i=0; i<this->size(); i++)
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
  if (this->size() == 1 && other.size()==1)
    {
      if (!(*this)[0].isComplex() && !other[0].isComplex())
        {
          return (*this)[0].sum(other[0], -1);
        }
      else return new ComplexExpr((*this)[0].real() - other[0].real(),
                                  (*this)[0].imag() - other[0].imag());
    }

  /* otherwise, create a list of the sums */
  Array<Expr> rtn(this->size());
  if (this->size() == other.size())
    {
      for (unsigned int i=0; i<this->size(); i++)
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
  if (this->size() == 1 && other.size()==1)
    {
      if (!(*this)[0].isComplex() && !other[0].isComplex())
        {
          return (*this)[0].multiply(other[0]);
        }
      else 
        {
          Expr me = (*this)[0];
          Expr you = other[0];
          if (Re(me).sameAs(Re(you)) && Im(me).sameAs(-Im(you)))
            {
              return Re(me)*Re(you) - Im(me)*Im(you);
            }
          else
            {
              return new ComplexExpr(Re(me)*Re(you) - Im(me)*Im(you),
                                     Re(me)*Im(you) + Im(me)*Re(you));
            }
        }
    }

  /* if the left operand is a scalar, multiply it through */
  if (this->size()==1)
    {
      Array<Expr> rtn(other.size());
      for (unsigned int i=0; i<other.size(); i++)
        {
          rtn[i] = (*this)[0] * other[i];
        }
      return new ListExpr(rtn);
    }

  /* if the right operand is a scalar, multiply it through */
  if (other.size()==1)
    {
      Array<Expr> rtn(this->size());
      for (unsigned int i=0; i<this->size(); i++)
        {
          rtn[i] = (*this)[i] * other[0];
        }
      return new ListExpr(rtn);
    }

  /* if both operands are flat vectors, take their dot product */
  if (this->size() == totalSize() && other.size()==other.totalSize() 
      && this->size() == other.size())
    {
      Expr rtn = new ZeroExpr();

      for (unsigned int i=0; i<this->size(); i++)
        {
          rtn = rtn + (*this)[i]*other[i];
        }
      return rtn;
    }

  /* if the left operand is a rectangular matrix and the 
   * right operand is a vector */
  unsigned int cols = (*this)[0].size();
  bool rectangular = true;
  for (unsigned int i=0; i<this->size(); i++)
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
  
  Array<Expr> rtn(this->size());
  for (unsigned int i=0; i<this->size(); i++)
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

  if (this->isComplex())
    {
      return new ComplexExpr(-real(), -imag());
    }

  /* if we are a scalar, process the unary minus here */
  if (this->size()==1)
    {
      const ConstantExpr* c = dynamic_cast<const ConstantExpr*>((*this)[0].ptr().get());
      const UnaryMinus* u = dynamic_cast<const UnaryMinus*>((*this)[0].ptr().get());
      /* if we are a constant, just modify the constant */
      if (c != 0)
        {
          if (c->value()==0.0)
            {
              return new ZeroExpr();
            }
          else
            {
              return new ConstantExpr(-1.0 * c->value());
            }
        }
      else if (u != 0) /* if we are already a unary minus, apply -(-x) --> x */
        {
          return u->arg();
        }
      else
        {
          RefCountPtr<ScalarExpr> sThis 
            = rcp_dynamic_cast<ScalarExpr>((*this)[0].ptr());
          TEST_FOR_EXCEPTION(sThis.get()==NULL, InternalError,
                             "Expr::operator-(): Operand " << (*this)[0].toString() 
                             << " is a non-scalar expression. All list structure "
                             "should have been handled before this point");
          return new UnaryMinus(sThis);
        }
    }

  /* otherwise, distribute the sign change over the list */
  Array<Expr> rtn(this->size());
  for (unsigned int i=0; i<this->size(); i++)
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

  /* If other is complex, transform to make the denominator real */
  if (other.isComplex())
    {
      Expr magSq = other.real()*other.real() + other.imag()*other.imag();
      return (*this) * other.conj() / magSq;
    }

  /* If I'm complex and the other is not, distribute division over re and im */
  if (isComplex() && !other.isComplex())
    {
      return new ComplexExpr(real()/other, imag()/other);
    }

  /* if we are a scalar, do simple scalar division */
  if (this->size()==1)
    {
      return (*this)[0].divide(other[0]);
    }

  /* otherwise, divide each element of the left by the right operand */
  Array<Expr> rtn(this->size());
  for (unsigned int i=0; i<this->size(); i++)
    {
      rtn[i] = (*this)[i] / other;
    }
  return new ListExpr(rtn);
}

const Expr& Expr::operator[](int i) const
{
  TEST_FOR_EXCEPTION(ptr().get()==NULL, InternalError,
                     "null this detected in Expr::operator[].");

  TEST_FOR_EXCEPTION(i<0 || i>=(int) this->size(), RuntimeError,
                     "Expr::operator[]() index i=" << i << " out of range [0, "
                     << this->size() << " in expr " << toString());

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

unsigned int Expr::size() const
{
  if (ptr().get()==0) return 0;

  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      return le->size();
    }
  return 1;
}

unsigned int Expr::totalSize() const
{
  if (ptr().get()==0) return 0;

  const ListExpr* le = dynamic_cast<const ListExpr*>(ptr().get());

  if (le != 0)
    {
      return le->totalSize();
    }
  return 1;
}

Expr Expr::real() const
{
  if (isComplex()) 
    {
      const ComplexExpr* cx = dynamic_cast<const ComplexExpr*>(ptr().get());
      return cx->real();
    }
  else if (size() > 1)
    {
      Array<Expr> rtn(size());
      for (unsigned int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i].real();
        }
      return new ListExpr(rtn);
    }
  else
    {
      return *this;
    }
}

Expr Expr::imag() const
{
  if (isComplex()) 
    {
      const ComplexExpr* cx = dynamic_cast<const ComplexExpr*>(ptr().get());
      return cx->imag();
    }
  else if (size() > 1)
    {
      Array<Expr> rtn(size());
      for (unsigned int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i].imag();
        }
      return new ListExpr(rtn);
    }
  else
    {
      return 0.0;
    }
  
}

Expr Expr::conj() const
{
  if (size()==1)
    {
      if (isComplex()) 
        {
          return new ComplexExpr(real(), -imag());
        }
      else return real();
    }
  else
    {
      Array<Expr> rtn(size());
      for (unsigned int i=0; i<size(); i++)
        {
          rtn[i] = (*this)[i].conj();
        }
      return new ListExpr(rtn);
    }
}

void Expr::setParameterValue(const double& value)
{
  Parameter* pe = dynamic_cast<Parameter*>(ptr().get());
  TEST_FOR_EXCEPTION(pe==0, RuntimeError, 
                     "Expr " << *this << " is not a Parameter expr, and "
                     "so setParameterValue() should not be called");
  pe->setValue(value);
}

double Expr::getParameterValue() const 
{
  const Parameter* pe = dynamic_cast<const Parameter*>(ptr().get());
  TEST_FOR_EXCEPTION(pe==0, RuntimeError, 
                     "Expr " << *this << " is not a Parameter expr, and "
                     "so getParameterValue() should not be called");
  return pe->value();
}

namespace SundanceCore
{
  using namespace SundanceUtils;

  Expr Complex(const Expr& re, const Expr& im)
  {
    TEST_FOR_EXCEPTION(re.size() != im.size(), RuntimeError,
                       "arguments mismatched in Complex(). Real part="
                       << re << ", imaginary part=" << im);

    TEST_FOR_EXCEPTION(re.isComplex() || im.isComplex(), RuntimeError,
                       "recursively defined complex number. Real part="
                       << re << ", imaginary part=" << im);

    if (re.totalSize() > 1U)
      {
        Array<Expr> rtn(re.size());
        for (unsigned int i=0; i<re.size(); i++)
          {
            rtn[i] = Complex(re[i], im[i]);
          }
        return new ListExpr(rtn);
      }

    const ZeroExpr* zr = dynamic_cast<const ZeroExpr*>(re[0].ptr().get());
    const ZeroExpr* zi = dynamic_cast<const ZeroExpr*>(im[0].ptr().get());

    if (zr == 0) /* nonzero real part */
      {
        if (zi==0) /* nonzero imag part */
          {
            return new ComplexExpr(re, im);
          }
        else /* zero imag part */
          {
            return re;
          }
      }
    else /* zero real part */
      {
        if (zi != 0) /* both are zero */
          {
            return Expr(0.0);
          }
        else /* pure imaginary */
          {
            return new ComplexExpr(0.0, im);
          }
      }
    return new ComplexExpr(re, im);
  }

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

  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d, const Expr& e, const Expr& f,
            const Expr& g)
  {
    return new ListExpr(tuple(a,b,c,d,e,f,g));
  }

  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d, const Expr& e, const Expr& f,
            const Expr& g, const Expr& h)
  {
    return new ListExpr(tuple(a,b,c,d,e,f,g,h));
  }
}

