/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiffOp.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

DiffOp::DiffOp(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg)
  : UnaryExpr(arg), mi_(op), myCoordDeriv_(), requiredFunctions_()
{
  myCoordDeriv_ = new CoordDeriv(mi_.firstOrderDirection());
}

ostream& DiffOp::toText(ostream& os, bool /* paren */) const 
{
	os << "D[" << arg().toString() << ", " << mi_.toString() << "]";
	return os;
}

ostream& DiffOp::toLatex(ostream& os, bool /* paren */) const 
{
	os << "D^{" << mi_.toString() << "}" << arg().toLatex();
	return os;
}

XMLObject DiffOp::toXML() const 
{
	XMLObject rtn("DiffOp");
	rtn.addAttribute("m", mi_.toString());
  rtn.addChild(arg().toXML());

	return rtn;
}



bool DiffOp::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  TimeMonitor t(nonzeroDerivCheckTimer());
  hasNonzeroDerivCalls()++;

  if (derivHasBeenCached(d))
    {
      nonzeroDerivCacheHits()++;
      return getCachedDerivNonzeroness(d);
    }

  TimeMonitor t2(uncachedNonzeroDerivCheckTimer());

  bool rtn = false;

  if (d.order()==0)
    {
      /* check whether there is an explicit dependence on our spatial var */
      MultipleDeriv me;
      me.put(myCoordDeriv_);
      if (evaluatableArg()->hasNonzeroDeriv(me))
        {
          rtn = true;
          addDerivToCache(d, rtn);

          return rtn;
        }
      /* See if we depend on any other functions */
      SundanceUtils::Set<Deriv> dependencies;
      evaluatableArg()->getRoughDependencies(dependencies);
      rtn = (dependencies.size() > 0);
    }
  else
    {
      /* check to see if the argument has a derivative wrt 
       * the set union {d, myCoordDeriv_} */
      MultipleDeriv dAndX = d;
      dAndX.put(myCoordDeriv_);
      if (evaluatableArg()->hasNonzeroDeriv(dAndX))
        {
          rtn = true;
          addDerivToCache(d, rtn);

          return rtn;
        }

      /* Check each term in the differentiated chain rule */
      SundanceUtils::Set<Deriv> argDeps;
      evaluatableArg()->getRoughDependencies(argDeps);

      Array<MultipleDeriv> leftOps;
      Array<MultipleDeriv> rightOps;

      d.productRulePermutations(leftOps, rightOps);

      SundanceUtils::Set<Deriv>::const_iterator iter;

      for (iter=argDeps.begin(); iter != argDeps.end(); iter++)
        {
          const Deriv& dep = *iter;
          
          /* Paranoiacally check that this dependency is a functional deriv */
          TEST_FOR_EXCEPTION(!dep.isFunctionalDeriv(), InternalError,
                             "DiffOp::hasNonzeroDeriv() detected "
                             "a non-functional derivative in its "
                             "dependency set");
                             
          const FunctionalDeriv* fDep = dep.funcDeriv();
          const MultiIndex& beta = fDep->multiIndex();
          int depFuncID = fDep->funcID();

          for (int j=0; j<leftOps.size(); j++)
            {
              if (rightOps[j].order() > 1) continue;

              if (rightOps[j].order()==1)
                {
                  Deriv r = *(rightOps[j].begin());
                  if (!r.isFunctionalDeriv()) continue;
                  const FunctionalDeriv* fr = r.funcDeriv();
                  const MultiIndex& gamma = fr->multiIndex();
                  int funcID = fr->funcID();
                  if (funcID != depFuncID || !(mi_+beta == gamma)) continue;
                }

              MultipleDeriv L = leftOps[j];
              L.put(dep);
              if (evaluatableArg()->hasNonzeroDeriv(L))
                {
                  rtn = true;
                  addDerivToCache(d, rtn);

                  return rtn;
                }
            }
        }
      
    }
  
  addDerivToCache(d, rtn);

  return rtn;
}



void DiffOp::getRoughDependencies(Set<Deriv>& funcs) const
{
  SundanceUtils::Set<Deriv> argDeps;
  evaluatableArg()->getRoughDependencies(argDeps);

  SundanceUtils::Set<Deriv>::const_iterator iter;

  for (iter=argDeps.begin(); iter != argDeps.end(); iter++)
    {
      const Deriv& d = *iter;
      if (d.isFunctionalDeriv())
        {
          funcs.put(d);
          funcs.put(d.funcDeriv()->derivWrtMultiIndex(mi_));
        }
      
    }
}

Array<DerivSet> DiffOp::derivsRequiredFromOperands(const DerivSet& derivs) const
{
  Tabs tabs;

  if (verbosity() > 1)
    {
      cerr << tabs << "Getting derivs required to evaluate diff op" << endl;
      cerr << tabs << "op = " << mi().toString() << endl;
      cerr << tabs << "arg = " << arg().toString() << endl;
    }

  DerivSet::const_iterator i;

  DerivSet nonzeros;
  MultipleDeriv d0;
  nonzeros.put(d0);

  for (i=derivs.begin(); i != derivs.end(); i++)
    {
      Tabs tabs0;

      const MultipleDeriv& d = *i;

      if (verbosity() > 0)
        {
          cerr << tabs0 << "funcs required to eval " << d << endl;
          cerr << tabs0 << "first trying mixed spatial-functional deriv" 
               << endl;
        }

      /* first check to see if (D_d D_{x_me})(arg) exists */
      MultipleDeriv dAndX = d;
      dAndX.put(myCoordDeriv_);
      if (evaluatableArg()->hasNonzeroDeriv(dAndX)) 
        {
          Tabs tabs1;
          if (verbosity() > 1)
            {
              cerr << tabs1 << "mixed deriv " << dAndX 
                   << " is nonzero " << endl;
            }
          nonzeros.put(dAndX);
        }
      else
        {
          Tabs tabs1;
          if (verbosity() > 1)
            {
              cerr << tabs1 << "mixed deriv " << dAndX << " is zero " << endl;
            }
        }

      /* Now try every term in the chain rule */
      SundanceUtils::Set<Deriv> argDeps;
      evaluatableArg()->getRoughDependencies(argDeps);

      Array<MultipleDeriv> leftOps;
      Array<MultipleDeriv> rightOps;

      d.productRulePermutations(leftOps, rightOps);

      SundanceUtils::Set<Deriv>::const_iterator j;

      for (j = argDeps.begin(); j != argDeps.end(); j++)
        {
          const Deriv& dep = *j;

          /* Paranoiacally check that this dependency is a functional deriv */
          TEST_FOR_EXCEPTION(!dep.isFunctionalDeriv(), InternalError,
                             "DiffOp::hasNonzeroDeriv() detected "
                             "a non-functional derivative in its "
                             "dependency set");

          const FunctionalDeriv* fDep = dep.funcDeriv();
          const MultiIndex& beta = fDep->multiIndex();
          int depFuncID = fDep->funcID();

          for (int k=0; k<leftOps.size(); k++)
            {
              if (rightOps[k].order() > 1) continue;

              if (rightOps[k].order()==1)
                {
                  Deriv r = *(rightOps[k].begin());
                  if (!r.isFunctionalDeriv()) continue;
                  const FunctionalDeriv* fr = r.funcDeriv();
                  const MultiIndex& gamma = fr->multiIndex();
                  int funcID = fr->funcID();
                  if (funcID != depFuncID || !(mi_+beta == gamma)) continue;
                }
              MultipleDeriv L = leftOps[k];
              L.put(dep);
              if (evaluatableArg()->hasNonzeroDeriv(L))
                {
                  if (rightOps[k].order()==0)
                    {
                      if (dep.isFunctionalDeriv()) 
                        {
                          requiredFunctions_[d].put(fDep->derivWrtMultiIndex(mi_));
                          if (verbosity() > 1)
                            {
                              Tabs tabs1;
                              cerr << tabs1 << "diff op will require " 
                                   << "evaluation of function " 
                                   << fDep->derivWrtMultiIndex(mi_).toString()
                                   << endl;
                            }
                        }
                    }
                  if (verbosity() > 1)
                    {
                      Tabs tabs1;
                      cerr << tabs1 << "found nonzero deriv " << L.toString()
                           << endl;
                    }
                  nonzeros.put(L);
                }
            }
        }
    }

  if (verbosity() > 1)
    {
      cerr << tabs << "operand's derivs are: " << endl;
      cerr << tabs << nonzeros << endl;
      cerr << tabs << "done getting required derivs of operand" << endl;
    }
  
  return tuple(nonzeros);
}


