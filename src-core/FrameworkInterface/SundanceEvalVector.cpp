/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvalVector.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


EvalVector::EvalVector()
  : LoadableVector(), 
    vectorVal_(),
    stringVal_(),
    constantVal_(0.0), 
    isConstant_(true),
    isZero_(true),
    isOne_(false)
{}

EvalVector::EvalVector(int n)
  : LoadableVector(), 
    vectorVal_(n),
    stringVal_(),
    constantVal_(0.0), 
    isConstant_(false),
    isZero_(false),
    isOne_(false)
{}

void EvalVector::setToZero() 
{
  constantVal_ = 0.0;
  isOne_ = false;
  isConstant_ = true;
  isZero_ = true;
}

void EvalVector::setToOne() 
{
  constantVal_ = 1.0;
  isOne_ = true;
  isConstant_ = true;
  isZero_ = false;
}

void EvalVector::setToConstantValue(const double& constantVal)
{
  constantVal_ = constantVal;
  isOne_ = false;
  isConstant_ = true;
  isZero_ = false;
}

void EvalVector::setToVectorValue()
{
  isOne_ = false;
  isConstant_ = false;
  isZero_ = false;
}

void EvalVector::setStringValue(const string& stringVal)
{
  stringVal_ = stringVal;
}

string EvalVector::getStringValue() const 
{
  if (isZero_) return "0.0";
  else if (isOne_) return "1.0";
  else if (isConstant_) return Teuchos::toString(constantVal_);
  else if (stringVal_.length()==0) 
    {
      if (verbosity() > VerbHigh)
        {
          return vectorVal_.toString();
        }
      else
        {
          return "{vector}";
        }
    }
  else return stringVal_;
}





void EvalVector::addScaled(const RefCountPtr<EvalVector>& other,
                           const double& scalar)
{
  Tabs tabs;
  if (verbosity() > VerbLow) 
    {
      cerr << tabs << "adding vectors " 
           << getStringValue() 
           << " and " << other->getStringValue() 
           << " with scale " << scalar << endl;
    }
  if (other->isZero() || scalar==0.0) 
    {
      if (verbosity() > VerbLow) cerr << tabs << "right=0" << endl;
      return;
    }

  

  if (isZero()) 
    {
      if (verbosity() > VerbLow) cerr << tabs << "left=0" << endl;
      if (other->isConstant())
        {
          setToConstantValue(other->getConstantValue() * scalar);
        }
      else 
        {
          resize(other->length());
          TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                             "mismatched vector lengths: me="
                             << length() << ", you=" << other->length());
          double* const x = start();
          const double* const y = other->start();
          if (scalar==1.0)
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] = y[i];
                }
              if (shadowOps()) setStringValue(other->getStringValue());
            }
          else if (scalar==-1.0)
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] = -y[i];
                }
              if (shadowOps()) setStringValue("-"+other->getStringValue());
            }
          else 
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] = scalar*y[i];
                }
              if (shadowOps()) setStringValue(Teuchos::toString(scalar)
                                              + "*" 
                                              + other->getStringValue());
            }
          setToVectorValue();
 
        }
    }
  else if (isConstant() && other->isConstant())
    {
      if (verbosity() > VerbLow) cerr << tabs << "both are constant" << endl;
      setToConstantValue(getConstantValue() 
                         + scalar*other->getConstantValue());
    }
  else if (isConstant())
    {
      if (verbosity() > VerbLow) cerr << tabs << "left is constant" << endl;
      resize(other->length());
      TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                         "mismatched vector lengths: me="
                         << length() << ", you=" << other->length());
      double c = getConstantValue();
      double* const x = start();
      const double* const y = other->start();
      if (scalar==1.0)
        {
          for (int i=0; i<length(); i++)
            {
              x[i] = c + y[i]; 
            }
          if (shadowOps()) 
            {
              setStringValue("(" + Teuchos::toString(getConstantValue())
                             + "+" + other->getStringValue() + ")");
            }
        }
      else if (scalar==-1.0)
        {
          for (int i=0; i<length(); i++)
            {
              x[i] = c - y[i]; 
            }
          if (shadowOps()) 
            {
              setStringValue("(" + Teuchos::toString(getConstantValue())
                             + "-" + other->getStringValue() + ")");
            }
        }
      else 
        {
          for (int i=0; i<length(); i++)
            {
              x[i] = c + scalar*y[i]; 
            }
          if (shadowOps()) 
            {
              setStringValue("(" + Teuchos::toString(getConstantValue())
                             + "+" + Teuchos::toString(scalar)
                             + "*" + other->getStringValue() + ")");
            }
        }
      setToVectorValue();
       
    }
  else if (other->isConstant())
    {
      if (verbosity() > VerbLow) cerr << tabs << "right is constant" << endl;
      double c = scalar*other->getConstantValue();

      double* const x = start();
      for (int i=0; i<length(); i++)
        {
          x[i] += c;
        }
      setToVectorValue();
      if (shadowOps()) setStringValue("(" + Teuchos::toString(getStringValue())
                                      + "+" + Teuchos::toString(c) + ")");
    }
  else
    {
      if (verbosity() > VerbLow) cerr << tabs << "both are non-constant" << endl;
      TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                         "mismatched vector lengths: me="
                         << length() << ", you=" << other->length());
      double* const x = start();
      const double* const y = other->start();
      if (scalar==1.0)
        {
          for (int i=0; i<length(); i++)
            {
              x[i] += y[i]; 
            }
          if (shadowOps()) setStringValue("(" + getStringValue() + "+" 
                                          + other->getStringValue() + ")");
        }
      else if (scalar==-1.0)
        {
          for (int i=0; i<length(); i++)
            {
              x[i] -= y[i]; 
            }
          if (shadowOps()) setStringValue("(" + getStringValue() + "-" 
                                          + other->getStringValue() +")");
        }
      else 
        {
          for (int i=0; i<length(); i++)
            {
              x[i] += scalar*y[i]; 
            }
          if (shadowOps()) setStringValue("(" + getStringValue() + "+" 
                                          + Teuchos::toString(scalar) 
                                          + "*" + other->getStringValue()
                                          +")");
        }
      setToVectorValue();
    }

  if (verbosity() > VerbLow) 
    {
      cerr << tabs << "result = " << getStringValue() << endl;
    }
}



void EvalVector::multiply(const RefCountPtr<EvalVector>& other)
{
  SUNDANCE_OUT(verbosity() > VerbLow, "multiplying " << getStringValue()
               << " and " << other->getStringValue());

  if (isZero() || other->isZero()) 
    {
      setToZero();
    }
  else if (other->isOne())
    {
      return;
    }
  else if (isConstant() && other->isConstant())
    {
      setToConstantValue(getConstantValue() * other->getConstantValue());
    }
  else if (isConstant())
    {
      resize(other->length());
      TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                         "mismatched vector lengths: me="
                         << length() << ", you=" << other->length());
      
      double c = getConstantValue();
      double* const x = start();
      const double* const y = other->start();
      for (int i=0; i<length(); i++)
        {
          x[i] = c*y[i];
        }
      setToVectorValue();
      if (shadowOps())
        {
          setStringValue(Teuchos::toString(getConstantValue())
                         + "*" + other->getStringValue());
        }
    }
  else if (other->isConstant())
    {
      double* const x = start();
      double c = other->getConstantValue();      
      for (int i=0; i<length(); i++)
        {
          x[i] *= c;
        }
      setToVectorValue();
      if (shadowOps())
        {
          setStringValue(Teuchos::toString(other->getConstantValue())
                         + "*" + getStringValue());
        }
    }
  else
    {
      TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                         "mismatched vector lengths: me="
                         << length() << ", you=" << other->length());                 
      double* const x = start();
      const double* const y = other->start();\
      for (int i=0; i<length(); i++)
        {
          x[i] *= y[i]; 
        }
      setToVectorValue();
      if (shadowOps())
        {
          setStringValue(getStringValue() + "*" + other->getStringValue());
        }
    }
}

void EvalVector::addProduct(const RefCountPtr<EvalVector>& a,
                            const RefCountPtr<EvalVector>& b)
{
  SUNDANCE_OUT(verbosity() > VerbLow, "adding product L=" << 
               a->getStringValue() << " R=" << b->getStringValue()
               << " to LHS " << getStringValue());

  if (a->isZero() || b->isZero()) 
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "one or both of the operands is zero, doing nothing");
      return;
    }

  if (a->isConstant())
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "first operand is constant, delegating to addScaled()");
      addScaled(b, a->getConstantValue());
      return;
    }
  
  if (b->isConstant())
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "first operand is constant, delegating to addScaled()");
      addScaled(a, b->getConstantValue());
      return;
    }

  resize(a->length());

  TEST_FOR_EXCEPTION(length()!=a->length(), InternalError,
                     "mismatched vector lengths: me="
                     << length() << ", a=" << a->length());
  TEST_FOR_EXCEPTION(length()!=b->length(), InternalError,
                     "mismatched vector lengths: me="
                     << length() << ", b=" << b->length());

  double* const x = start();
  const double* const y = a->start();
  const double* const z = b->start();
      
  if (isZero())
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "initial value is zero, adding vector dot-slash product");
      for (int i=0; i<length(); i++)
        {
          x[i] = y[i]*z[i];
        }
      if (shadowOps()) setStringValue("(" + a->getStringValue() 
                                      + "*" + b->getStringValue() 
                                      + ")");
    }
  else if (isConstant())
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "initial value is constant, adding to vector dot-slash product");
      double c = getConstantValue();
      for (int i=0; i<length(); i++)
        {
          x[i] = c + y[i]*z[i];
        }
      if (shadowOps()) setStringValue("(" + Teuchos::toString(c) 
                                      + "+(" + a->getStringValue() 
                                      + "*" + b->getStringValue() 
                                      + "))");
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "initial value is a vector, adding to vector dot-slash product");
      for (int i=0; i<length(); i++)
        {
          x[i] += y[i]*z[i];
        }
      if (shadowOps()) setStringValue("(" + getStringValue() + "+(" 
                                      + a->getStringValue() 
                                      + "*" + b->getStringValue() 
                                      + "))");
    }
  setToVectorValue();
}

void EvalVector::sqrt()
{
  if (isConstant()) 
    {
      setToConstantValue(::sqrt(getConstantValue()));
      return;
    }
  else
    {
      double* const x = start();
      
      for (int i=0; i<length(); i++)
        {
          x[i] = ::sqrt(x[i]);
        }
      
      if (shadowOps()) setStringValue("sqrt(" + getStringValue() + ")");
    }
}

void EvalVector::unaryMinus()
{
  if (isConstant()) 
    {
      setToConstantValue(-1.0*getConstantValue());
    }
  else
    {
      double* const x = start();
          
      for (int i=0; i<length(); i++)
        {
          x[i] *= -1.0;
        }

      if (shadowOps()) setStringValue("(-" + getStringValue() + ")");
    }
}

void EvalVector::applyUnaryFunction(const UnaryFunctor* func,
                                    Array<RefCountPtr<EvalVector> >& funcDerivs) const 
{
  for (int i=0; i<funcDerivs.size(); i++)
    {
      funcDerivs[i] = rcp(new EvalVector());
      funcDerivs[i]->setToVectorValue();
    }
  
  if (isConstant()) 
    {
      double xVal = getConstantValue();
      const double* const x = &xVal;
      double f;
      double df;
      if (funcDerivs.size()==1)
        {
          func->eval(x, 1, &f);
          funcDerivs[0]->setToConstantValue(f);
        }
      else
        {
          func->eval(x, 1, &f, &df);
          funcDerivs[0]->setToConstantValue(f);
          funcDerivs[1]->setToConstantValue(f);
        }
    }
  else
    {
      const double* const x = start();
      if (funcDerivs.size()==1)
        {
          funcDerivs[0]->resize(length());
          double* f = funcDerivs[0]->start();
          func->eval(x, length(), f);
        }
      else
        {
          funcDerivs[0]->resize(length());
          funcDerivs[1]->resize(length());
          double* f = funcDerivs[0]->start();
          double* df = funcDerivs[1]->start();
          func->eval(x, length(), f, df);
        }
      
    }
      
  if (shadowOps())
    { 
      for (int i=0; i<funcDerivs.size(); i++)
        {
          if (i==0) 
            {
              funcDerivs[i]->setStringValue(func->name() 
                                            + "(" + getStringValue() + ")");
            }
          else
            {
              string primes;
              for (int j=1; j<=i; j++) primes += "'";
              funcDerivs[i]->setStringValue(func->name() + primes
                                            + "(" + getStringValue() + ")");
            }
        }
    }
}

void EvalVector::copy(const RefCountPtr<EvalVector>& other)
{
  if (other->isZero()) 
    {
      setToZero();
    }
  else if (other->isOne())
    {
      setToOne();
    }
  else if (other->isConstant())
    {
      setToConstantValue(other->getConstantValue());
    }
  else 
    {
      resize(other->length());


      TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                         "mismatched vector lengths in EvalVector::copy()");

      double* const x = start();
      const double* const y = other->start();

      for (int i=0; i<length(); i++)
        {
          x[i] = y[i];
        }
      setToVectorValue();
      if (shadowOps()) setStringValue(other->getStringValue());
    }
}

void EvalVector::print(ostream& os) const 
{
  os << "[";
  if (isConstant_) os << constantVal_;
  else
    {
      if (shadowOps()) os << stringVal_;
      else 
        {
          if (verbosity() > VerbHigh)
            {
              os << vectorVal_;
            }
          else
            {
              os << "{vector}" << endl;
            }
        }
    }
  os << "]";
}
