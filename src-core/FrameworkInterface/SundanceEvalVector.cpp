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
    isOne_(false),
    numerical_(true)
{}

EvalVector::EvalVector(int n)
  : LoadableVector(), 
    vectorVal_(n),
    stringVal_(),
    constantVal_(0.0), 
    isConstant_(false),
    isZero_(false),
    isOne_(false),
    numerical_(true)
{}

void EvalVector::setToZero() 
{
  constantVal_ = 0.0;
  isOne_ = false;
  isConstant_ = true;
  isZero_ = true;
  numerical_ = true;
}

void EvalVector::setToOne() 
{
  constantVal_ = 1.0;
  isOne_ = true;
  isConstant_ = true;
  isZero_ = false;
  numerical_ = true;
}

void EvalVector::setToConstantValue(const double& constantVal)
{
  constantVal_ = constantVal;
  isOne_ = false;
  isConstant_ = true;
  isZero_ = false;
  numerical_ = true;
}

void EvalVector::setToVectorValue()
{
  isOne_ = false;
  isConstant_ = false;
  isZero_ = false;
  numerical_ = true;
}

void EvalVector::setToStringValue(const string& stringVal)
{
  isOne_ = false;
  isConstant_ = false;
  isZero_ = false;
  stringVal_ = stringVal;
  numerical_ = false;
}

string EvalVector::getStringValue() const 
{
  if (isZero_) return "0.0";
  else if (isOne_) return "1.0";
  else if (isConstant_) return Teuchos::toString(constantVal_);
  else if (stringVal_.length()==0) return vectorVal_.toString();
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
          if (numerical())
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
                }
              else if (scalar==-1.0)
                {
                  for (int i=0; i<length(); i++)
                    {
                      x[i] = -y[i];
                    }
                }
              else 
                {
                  for (int i=0; i<length(); i++)
                    {
                      x[i] = scalar*y[i];
                    }
                }
              setToVectorValue();
            }
          else
            {
              if (scalar==1.0)
                {
                  setToStringValue(other->getStringValue());
                }
              else if (scalar==-1.0)
                {
                  setToStringValue("-" + other->getStringValue());
                }
              else 
                {
                  setToStringValue(Teuchos::toString(scalar)
                                 + "*" + other->getStringValue());
                }
            }
 
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
      if (numerical())
        {
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
            }
          else if (scalar==-1.0)
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] = c - y[i]; 
                }
            }
          else 
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] = c + scalar*y[i]; 
                }
            }
          setToVectorValue();
        }
      else /* phony calculation with strings */
        {
          if (scalar==1.0)
            {
              setToStringValue("(" + Teuchos::toString(getConstantValue())
                             + "+" + other->getStringValue() + ")");
            }
          else if (scalar==-1.0)
            {
              setToStringValue("(" + Teuchos::toString(getConstantValue())
                             + "-" + other->getStringValue() + ")");
            }
          else
            {
              setToStringValue("(" + Teuchos::toString(getConstantValue())
                             + "+" + Teuchos::toString(scalar)
                             + "*" + other->getStringValue() + ")");
            }
        }
    }
  else if (other->isConstant())
    {
      if (verbosity() > VerbLow) cerr << tabs << "right is constant" << endl;
      double c = scalar*other->getConstantValue();

      if (numerical())
        {
          double* const x = start();
          for (int i=0; i<length(); i++)
            {
              x[i] += c;
            }
          setToVectorValue();
        }
      else
        {
          setToStringValue("(" + Teuchos::toString(getStringValue())
                         + "+" + Teuchos::toString(c) + ")");
        }
    }
  else
    {
      if (verbosity() > VerbLow) cerr << tabs << "both are non-constant" << endl;
      if (numerical())
        {
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
            }
          else if (scalar==-1.0)
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] -= y[i]; 
                }
            }
          else 
            {
              for (int i=0; i<length(); i++)
                {
                  x[i] += scalar*y[i]; 
                }
            }
          setToVectorValue();
        }
      else
        {
          if (scalar==1.0)
            {
              setToStringValue("(" + getStringValue() + "+" 
                             + other->getStringValue() + ")");
            }
          else if (scalar==-1.0)
            {
              setToStringValue("(" + getStringValue() + "-" 
                             + other->getStringValue() +")");
            }
          else 
            {
              setToStringValue("(" + getStringValue() + "+" 
                             + Teuchos::toString(scalar) 
                             + "*" + other->getStringValue()+")");
            }
        }
    }

  if (verbosity() > VerbLow) 
    {
      cerr << tabs << "result = " << getStringValue() << endl;
    }
}



void EvalVector::multiply(const RefCountPtr<EvalVector>& other)
{
  cerr << "MULTIPLYING ***************************************" << endl;
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
      if (numerical())
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
        }
      else
        {
          setToStringValue(Teuchos::toString(getConstantValue())
                           + "*" + other->getStringValue());
        }
    }
  else if (other->isConstant())
    {
      if (numerical())
        {
          double* const x = start();
          double c = other->getConstantValue();      
          for (int i=0; i<length(); i++)
            {
              x[i] *= c;
            }
          setToVectorValue();
        }
      else
        {
          setToStringValue(Teuchos::toString(other->getConstantValue())
                           + "*" + getStringValue());
        }
    }
  else
    {
      if (numerical())
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
        }
      else
        {
          setToStringValue(getStringValue() + "*" + other->getStringValue());
        }
    }
}

void EvalVector::addProduct(const RefCountPtr<EvalVector>& a,
                            const RefCountPtr<EvalVector>& b)
{
  SUNDANCE_OUT(verbosity() > VerbLow, "adding product L=" << 
               a->getStringValue() << " R=" << b->getStringValue()
               << " to LHS " << getStringValue());

  if (a->isZero() || b->isZero()) return;

  if (a->isConstant())
    {
      addScaled(b, a->getConstantValue());
      return;
    }
  
  if (b->isConstant())
    {
      addScaled(a, b->getConstantValue());
      return;
    }

  if (numerical())
    {
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
          for (int i=0; i<length(); i++)
            {
              x[i] = y[i]*z[i];
            }
          setToVectorValue();
        }
      else if (isConstant())
        {
          double c = getConstantValue();
          for (int i=0; i<length(); i++)
            {
              x[i] = c + y[i]*z[i];
            }
          setToVectorValue();
        }
      else
        {
          for (int i=0; i<length(); i++)
            {
              x[i] += y[i]*z[i];
            }
        }
    }
  else
    {
      if (isZero())
        {
          setToStringValue("(" + a->getStringValue() 
                           + "*" + b->getStringValue() + ")");
        }
      else if (isConstant())
        {
          double c = getConstantValue();
          setToStringValue("(" + Teuchos::toString(c) 
                           + "+(" + a->getStringValue() 
                           + "*" + b->getStringValue() + "))");
        }
      else
        {
          setToStringValue("(" + getStringValue() + "+(" 
                           + a->getStringValue() 
                           + "*" + b->getStringValue() + "))");
        }
    }
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
      if (numerical())
        {
          double* const x = start();
          
          for (int i=0; i<length(); i++)
            {
              x[i] = ::sqrt(x[i]);
            }
        }
      else
        {
          setToStringValue("sqrt(" + getStringValue() + ")");
        }
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
      if (numerical())
        {
          double* const x = start();
          
          for (int i=0; i<length(); i++)
            {
              x[i] *= -1.0;
            }
        }
      else
        {
          setToStringValue("(-" + getStringValue() + ")");
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
      if (numerical())
        {
          resize(other->length());

          TEST_FOR_EXCEPTION(length()!=other->length(), InternalError,
                             "mismatched vector lengths");

          double* const x = start();
          const double* const y = other->start();

          for (int i=0; i<length(); i++)
            {
              x[i] = y[i];
            }
          setToVectorValue();
        }
      else
        {
          setToStringValue(other->getStringValue());
        }
    }
}

void EvalVector::print(ostream& os) const 
{
  if (numerical()) 
    {
      os << "numerical[";
      if (isConstant_) os << constantVal_;
      else os << vectorVal_;
      os << "]";
    }
  else
    {
      os << "string[" << stringVal_ << "]" << endl;
    }
}
