/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace TSFExtended;
using namespace Internal;
using namespace Internal;

static Time& derivSetLookupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("deriv set lookup"); 
  return *rtn;
}

EvaluatableExpr::EvaluatableExpr()
	: ScalarExpr(), 
    regionToDerivSetIndexMap_(),
    derivSetToDerivSetIndexMap_(),
    derivSets_(),
    sparsityPatterns_(),
    evaluators_(),
    derivNonzeronessCache_(),
    evalRefCount_(0),
    resultCacheIsValid_(false),
    resultCache_(),
    currentDerivSuperset_()
{}

bool EvaluatableExpr::hasWorkspace() const
{
  TEST_FOR_EXCEPTION(refCount() <= 0, InternalError,
                     "bad ref count in EvaluatableExpr::hasWorkspace()");
  refCount()--;
  if (refCount()==0) return true;
  return false;
}





DerivSet EvaluatableExpr::identifyNonzeroDerivs() const
{
  DerivSet nonzeroDerivs;

  /* the zeroth deriv is always nonzero */
  nonzeroDerivs.put(MultipleDeriv());

  /* Look for first derivs */
  SundanceUtils::Set<Deriv> d;
  getRoughDependencies(d);



  Array<Deriv> d1;
  
  for (Set<Deriv>::const_iterator i=d.begin(); i != d.end(); i++)
    {
      d1.append(*i);
      MultipleDeriv m;
      m.put(*i);
      if (hasNonzeroDeriv(m)) 
        {
          nonzeroDerivs.put(m);
        }
    }
  /* assemble second derivs*/
  for (int i=0; i<d1.size(); i++)
    {
      for (int j=0; j<=i; j++)
        {
          MultipleDeriv m;
          m.put(d1[i]);
          m.put(d1[j]);
          if (hasNonzeroDeriv(m)) 
            {
              nonzeroDerivs.put(m);
            }
        }
    }

  return nonzeroDerivs;
  
}

bool EvaluatableExpr::checkForKnownRegion(const EvalContext& region,
                                          bool& derivSetIsKnown) const
{

  Tabs tabs;
  if (verbosity() > 1) 
    {
      cerr << tabs << "checking for region " << region.toString()
           << " in expr " 
           << toString() << endl;
    }
  if (regionToDerivSetIndexMap_.containsKey(region))
    {
      if (verbosity() > 1)
        {
          Tabs tabs1;
          cerr << tabs1 << "region " << region.toString() 
               << " already known "
               << endl;
        }
      if (derivSetToDerivSetIndexMap_.containsKey(currentDerivSuperset()))
        {
          derivSetIsKnown = true;
        }
      else
        {
          derivSetIsKnown = false;
        }
      return true;
    }
  else
    {
      Tabs tabs1;
      if (verbosity() > 1) 
      {
        cerr << tabs1 << "region " << region.toString() << " is new "
           << endl;
        cerr << tabs1 << "binding it to deriv set " << currentDerivSuperset()
             << endl;
      }
      /* if the region is new but the deriv set is already known,
       * return true */
      if (derivSetToDerivSetIndexMap_.containsKey(currentDerivSuperset()))
        {
          Tabs tabs2;

          int i = derivSetToDerivSetIndexMap_.get(currentDerivSuperset());
          if (verbosity() > 1)
            {
              cerr << tabs2 << "region is unknown, but deriv set " 
                   << currentDerivSuperset()
                   << " is already known, with index " << i << endl;
            }

          derivSetIsKnown = true;
        }
      else
        {
          /* If both the region and deriv set are new, return false,
           * which will tell the client that createNewDerivSet() needs
           * to be called */
          Tabs tabs2;
          if (verbosity() > 1)
            {
              cerr << tabs2 << "deriv set " << currentDerivSuperset()
                   << " is new " << endl;
            }
          derivSetIsKnown = false;
        }
      return false;
    }
}


int EvaluatableExpr::registerRegion(const EvalContext& region,
                                    bool derivSetIsKnown,
                                    const DerivSet& derivs,
                                    const EvaluatorFactory* factory,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;



  if (verbosity() > 1) 
    {
      cerr << tabs << "expr " << toString() << " registering deriv set "
           << derivs.toString() << endl;
    }

  /* paranoid check: make sure this region hasn't already been processed. */
  TEST_FOR_EXCEPTION(derivSetIsKnown 
                     && regionToDerivSetIndexMap_.containsKey(region),
                     InternalError,
                     "region " << region.toString() << " already registered "
                     "with expr " << toString() << 
                     " but EvaluatableExpr::registerRegion() was called "
                     "nonetheless");
  
  /* paranoid check: make sure this deriv set hasn't already been added */
  TEST_FOR_EXCEPTION(!derivSetIsKnown 
                     && derivSetToDerivSetIndexMap_.containsKey(derivs),
                     InternalError,
                     "derivSet " << derivs << " already registered with expr "
                     << toString() << 
                     " but EvaluatableExpr::registerRegion() was called anyway");

  int rtn;
  if (derivSetIsKnown)
    {
      /* if the deriv set is known, get its index */
      rtn = derivSetToDerivSetIndexMap_.get(derivs);
    }
  else
    {
      /* if this is a new deriv set, then
       * the index for this deriv set is the next index 
       * in the sparsity array */
      rtn = sparsityPatterns_.size();
    }

  if (verbosity() > 1) 
    {
      Tabs tabs1;
      cerr << tabs1 << "deriv set index will be " << rtn << endl;
    }
  
  /* If the deriv set is new, 
   * create the sparsity pattern for this deriv set */
  if (!derivSetIsKnown)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, tabs << "creating sparsity pattern");
      sparsityPatterns_.append(rcp(new SparsityPattern(derivs, this, regardFuncsAsConstant)));
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, tabs << "reusing sparsity pattern");
    }

  SUNDANCE_OUT(verbosity() > VerbMedium,
               tabs << "found sparsity pattern:\n"
               << *(sparsityPatterns_[rtn]));


  /** If the deriv set is new, 
   * create the evaluator for this deriv set */
  if (!derivSetIsKnown)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, tabs << "creating evaluator");
      evaluators_.append(rcp(factory->createEvaluator(this, rtn)));
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, tabs << "reusing evaluator");
    }

  /* set up mapping from deriv set to deriv set index */  
  if (!derivSetIsKnown)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tabs << "registering deriv set at DSI=" << rtn);
      derivSetToDerivSetIndexMap_.put(derivs, rtn);
      derivSets_.append(derivs);
    }

  regionToDerivSetIndexMap_.put(region, rtn);

  return rtn;
}



void EvaluatableExpr::flushResultCache() const 
{
  if (verbosity() > 2) 
    {
      cerr << "flushing result cache for " << toString() << endl;
    }
  resultCacheIsValid_=false;
}

void EvaluatableExpr::evaluate(const EvalManager& mgr,
                               RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor t(evalTimer());

  Tabs tab;

  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "evaluating " << toString());



  if (!resultCacheIsValid_)
    {
      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   "computing values for " << toString());
      int derivSetIndex = getDerivSetIndex(mgr.getRegion());
      evaluator(derivSetIndex)->eval(mgr, resultCache_);
      resultCacheIsValid_ = true;
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbMedium, 
                   "reusing values for " << toString());
    }
  results = resultCache_;
}

void EvaluatableExpr::findDerivSuperset(const DerivSet& derivs) const 
{
  Tabs tabs;
  if (verbosity() > 1)
    {
      Tabs tabs1;
      cerr << tabs << "finding deriv superset for expr: " << toString()
           << endl;
      cerr << tabs1 << "required derivs are " << derivs.toString() << endl;
      cerr << tabs1 << "current derivs superset is " 
           << currentDerivSuperset().toString() << endl;
    }
  currentDerivSuperset().merge(derivs);

  if (verbosity() > 1)
    {
      Tabs tabs1;
      cerr << tabs1 << "merged superset is " 
           << currentDerivSuperset().toString() << endl;
    }
}

int EvaluatableExpr::getDerivSetIndex(const EvalContext& region) const
{
  TimeMonitor timer(derivSetLookupTimer());
  Tabs tabs;
  if (verbosity() > 1)
    {
      cerr << tabs << "finding deriv set index for region: " 
           << region.toString()
           << endl;
      cerr << tabs << "map is " << regionToDerivSetIndexMap_ << endl;
    }

  TEST_FOR_EXCEPTION(!regionToDerivSetIndexMap_.containsKey(region),
                     InternalError,
                     "EvaluatableExpr::getDerivSetIndex: region " 
                     << region.toString() << " not found in map " 
                     << regionToDerivSetIndexMap_);

  int rtn = regionToDerivSetIndexMap_.get(region);

  if (verbosity() > 1)
    {
      cerr << tabs << "found DSI = " << rtn << endl;
    }

  return rtn;
  
}

const EvaluatableExpr* EvaluatableExpr::getEvalExpr(const Expr& expr)
{
  const EvaluatableExpr* rtn 
    = dynamic_cast<const EvaluatableExpr*>(expr[0].ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError,
                     "cast of " << expr 
                     << " failed in EvaluatableExpr::getEvalExpr()");
  TEST_FOR_EXCEPTION(expr.size() != 1, InternalError,
                     "non-scalar expression " << expr
                     << " in EvaluatableExpr::getEvalExpr()");

  return rtn;
}
