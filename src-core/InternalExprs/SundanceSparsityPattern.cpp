/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSparsityPattern.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceTabs.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;
using namespace TSFExtended;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SparsityPattern::SparsityPattern(const DerivSet& derivs,
                                 const EvaluatableExpr* expr,
                                 bool regardFuncsAsConstant)
  : derivToIndexMap_(),
    derivs_(),
    states_(),
    isFirstOrderSpatialDeriv_(),
    spatialDerivDir_()
{
  verbosity() = classVerbosity();
  DerivSet::const_iterator i;

  //  cerr << "computing sparsity pattern for " << expr->toString() << endl;

  SundanceUtils::Set<Deriv> singleDerivs;

  for (i=derivs.begin(); i != derivs.end(); i++)
    {
      const MultipleDeriv& d = *i;
      MultiSet<Deriv>::const_iterator j;
      for (j=d.begin(); j != d.end(); j++)
        {
          const Deriv& sd = *j;
          if (regardFuncsAsConstant)
            {
              if (!sd.isUnknownFunction()) singleDerivs.put(sd);
            }
          else
            {
              singleDerivs.put(sd);
            }
        }
    }

  if (!regardFuncsAsConstant)
    {
      DerivSet secondOrderDerivs = expr->identifyNonzeroDerivs();
      for (i=secondOrderDerivs.begin(); i != secondOrderDerivs.end(); i++)
        {
          const MultipleDeriv& d = *i;
          MultiSet<Deriv>::const_iterator j;
          for (j=d.begin(); j != d.end(); j++)
            {
              const Deriv& sd = *j;
              singleDerivs.put(sd);
            }
        }
    }

  for (int dir=0; dir<CoordDeriv::maxDim(); dir++)
    {
      singleDerivs.put(new CoordDeriv(dir));
    }

  for (i=derivs.begin(); i != derivs.end(); i++)
    {
      if (verbosity() > VerbHigh)
        {
          cerr << "============ looking for coeff of " << *i 
               << " in " << expr->toString() << endl;
        }
      const MultipleDeriv& d = *i;
      int n = derivs_.size();
      derivs_.append(d);
      derivToIndexMap_.put(d, n);
      if (!expr->hasNonzeroDeriv(d)) 
        {
          if (verbosity() > VerbHigh)
            {
              cerr << "deriv wrt " << d << " is zero" << endl;
            }
          states_.append(ZeroDeriv);
        }
      else  
        {
          if (expr->isConstant())
            {
              if (verbosity() > VerbHigh)
                {
                  cerr << "deriv " << d 
                       << " has identically constant coeff" << endl;
                }
              states_.append(ConstantDeriv);
            }
          else
            {
              /* to test whether a derivative is a constant, see
               * if its derivative wrt all of the single derivatives
               * is zero */
              SundanceUtils::Set<Deriv>::const_iterator j;
              bool isConstant = true;
              MultipleDeriv dTry;
              for (j = singleDerivs.begin(); j != singleDerivs.end(); j++)
                {
                  dTry = d;
                  dTry.put(*j);
                  if (verbosity() > VerbHigh)
                    {
                      cerr << "testing against " << dTry << endl;
                    }
                  if (expr->hasNonzeroDeriv(dTry)) 
                    {
                      if (verbosity() > VerbHigh)
                        {
                          cerr << "has nonzero deriv wrt" << dTry << endl;
                        }
                      isConstant = false;
                      break;
                    }
                  else
                    {
                      if (verbosity() > VerbHigh)
                        {
                          cerr << "has zero deriv wrt" << dTry << endl;
                        }
                    }
                }
              if (isConstant) 
                {
                  if (verbosity() > VerbHigh)
                    {
                      cerr << "deriv wrt " << d << " is constant" << endl;
                    }
                  states_.append(ConstantDeriv);
                }
              else 
                {
                  if (verbosity() > VerbHigh)
                    {
                      cerr << "deriv wrt " << d << " is non-constant" << endl;
                    }
                  states_.append(VectorDeriv);
                }
            }
        }
      if (d.order()==1 && d.begin()->isCoordDeriv())
        {
          isFirstOrderSpatialDeriv_.append(true);
          spatialDerivDir_.append(d.begin()->coordDeriv()->dir());
        }
      else
        {
          isFirstOrderSpatialDeriv_.append(false);
          spatialDerivDir_.append(-1);
        }
    }
}

bool SparsityPattern::containsDeriv(const MultipleDeriv& d) const
{
  return derivToIndexMap_.containsKey(d);
}

int SparsityPattern::getIndex(const MultipleDeriv& d) const
{
  return derivToIndexMap_.get(d);
}



void SparsityPattern::print(ostream& os) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  os << tabs << "SparsityPattern" << endl;
  for (int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\tderiv=\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs_[i].toString() << "\tstate=\t" ;
      switch(states_[i])
        {
        case ZeroDeriv:
          os  << "Zero" << endl;
          break;
        case ConstantDeriv:
          os << "Constant" << endl;
          break;
        case VectorDeriv:
          os << "Vector" << endl;
          break;
        }
    }
}

string SparsityPattern::toString() const 
{
	ostringstream ss;
	print(ss);
	return ss.str();
}
