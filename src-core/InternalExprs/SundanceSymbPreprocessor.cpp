/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace Internal;
using namespace Internal;

static Time& findNonzerosTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("identifying nonzero functional derivs "); 
  return *rtn;
}

DerivSet SymbPreprocessor::setupExpr(const Expr& expr, 
                                     const EvalContext& region, 
                                     const EvaluatorFactory* factory)
{
  TimeMonitor t(preprocTimer());

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  DerivSet derivs;

  derivs.put(MultipleDeriv());

  e->resetDerivSuperset();

  e->findDerivSuperset(derivs);

  bool regardFuncsAsConstant = false;
  e->setupEval(region, factory, regardFuncsAsConstant);

  return derivs;
}

DerivSet SymbPreprocessor::setupExpr(const Expr& expr, 
                                     const Expr& tests,
                                     const Expr& unks,
                                     const Expr& u0, 
                                     const EvalContext& region, 
                                     const EvaluatorFactory* factory,
                                     int maxDiffOrder)
{
  TimeMonitor t(preprocTimer());


  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  bool missingTestFunction = !e->allTermsHaveTestFunctions();
  TEST_FOR_EXCEPTION(missingTestFunction, RuntimeError,
                     "Weak form " << expr << " contains at least one "
                     "term without a test function");

  bool u0IsZero = false;
  DerivSet derivs = identifyNonzeroDerivs(expr, tests, unks, u0, maxDiffOrder,
                                          u0IsZero);


  e->resetDerivSuperset();

  e->findDerivSuperset(derivs);

  e->setupEval(region, factory, u0IsZero);

  return derivs;
}

DerivSet SymbPreprocessor::identifyNonzeroDerivs(const Expr& expr,
                                                 const Expr& tests,
                                                 const Expr& unks,
                                                 const Expr& evalPts,
                                                 int maxDiffOrder,
                                                 bool& u0IsZero)
{
  TimeMonitor t(findNonzerosTimer());


  /* first we have to make sure the expression is evaluatable */
  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(e==0, InternalError, 
                     "SymbPreprocessor::identifyNonzeroDerivs "
                     "given non-evaluatable expr" << expr.toString());

  /* make flat lists of tests and unknowns */
  Expr v = tests.flatten();
  Expr u = unks.flatten();
  Expr u0 = evalPts.flatten();

  SundanceUtils::Set<int> testID;
  SundanceUtils::Set<int> unkID;

  /* check the test functions for redundancies and non-test funcs */
  for (int i=0; i<v.size(); i++) 
    {
      const TestFuncElement* vPtr 
        = dynamic_cast<const TestFuncElement*>(v[i].ptr().get());
      TEST_FOR_EXCEPTION(vPtr==0, RuntimeError, "list of purported test funcs "
                         "contains a non-test function " << v[i].toString());
      int fid = vPtr->funcID();
      TEST_FOR_EXCEPTION(testID.contains(fid), RuntimeError,
                         "duplicate test function in list " 
                         << tests.toString());
      testID.put(fid);
      vPtr->substituteZero();
    }

  /* check the unk functions for redundancies and non-unk funcs */
  for (int i=0; i<u.size(); i++) 
    {
      const UnknownFuncElement* uPtr 
        = dynamic_cast<const UnknownFuncElement*>(u[i].ptr().get());
      TEST_FOR_EXCEPTION(uPtr==0, RuntimeError, 
                         "list of purported unknown funcs "
                         "contains a non-unknown function " 
                         << u[i].toString());
      int fid = uPtr->funcID();
      TEST_FOR_EXCEPTION(unkID.contains(fid), RuntimeError,
                         "duplicate unknown function in list " 
                         << unks.toString());
      unkID.put(fid);
      RefCountPtr<DiscreteFuncElement> u0Ptr 
        = rcp_dynamic_cast<DiscreteFuncElement>(u0[i].ptr());
      RefCountPtr<ZeroExpr> u0ZeroPtr 
        = rcp_dynamic_cast<ZeroExpr>(u0[i].ptr());
      TEST_FOR_EXCEPTION(u0Ptr.get()==NULL && u0ZeroPtr.get()==NULL, 
                         RuntimeError,
                         "evaluation point " << u0[i].toString() 
                         << " is neither a discrete function nor a zero expr");
      if (u0Ptr.get()==NULL)
        {
          u0IsZero = true;
          uPtr->substituteZero();
        }
      else
        {
          uPtr->substituteFunction(u0Ptr);
        }
    }

  

  DerivSet nonzeroDerivs;

  /* the zeroth deriv should be computed only if the max diff order is zero */
  if (maxDiffOrder==0) 
    {
      nonzeroDerivs.put(MultipleDeriv());
      return nonzeroDerivs;
    }
  

  /* Find any functions that might be in the expr */
  SundanceUtils::Set<Deriv> d;
  e->getRoughDependencies(d);

  SundanceUtils::Set<Deriv> testDerivs;
  SundanceUtils::Set<Deriv> unkDerivs;

  for (Set<Deriv>::const_iterator i=d.begin(); i != d.end(); i++)
    {
      const Deriv& di = *i;
      if (!di.isFunctionalDeriv()) continue;
      const FunctionalDeriv* fd = di.funcDeriv();
      int fid = fd->funcID();
      if (testID.contains(fid)) 
        {
          testDerivs.put(di);
          MultipleDeriv m1;
          m1.put(di);
          if (e->hasNonzeroDeriv(m1)) nonzeroDerivs.put(m1);
        }
      else if (unkID.contains(fid))
        {
          unkDerivs.put(di);
        }
      else
        {
          TEST_FOR_EXCEPTION(true, RuntimeError,
                             "Deriv " << di.toString() 
                             << " does not appear in either the test function "
                             "list " << v.toString() 
                             << " or the unknown function list " 
                             << u.toString());
        }
    }

  if (maxDiffOrder==2)
    {
      for (Set<Deriv>::const_iterator i=testDerivs.begin(); i != testDerivs.end(); i++)
        {
          const Deriv& dTest = *i;
          for (Set<Deriv>::const_iterator j=unkDerivs.begin(); j != unkDerivs.end(); j++)
            {
              const Deriv& dUnk = *j;
              
              MultipleDeriv m;
              m.put(dTest);
              m.put(dUnk);
              if (e->hasNonzeroDeriv(m))
                {
                  nonzeroDerivs.put(m);
                }
            }
        }
    }

  /* Sanity check: ensure there are no expressions that are nonlinear in a 
   * test function. */
  DerivSet badDerivs;
  for (Set<Deriv>::const_iterator i=testDerivs.begin(); i != testDerivs.end(); i++)
    {
      const Deriv& dTest1 = *i;
      for (Set<Deriv>::const_iterator j=testDerivs.begin(); j != testDerivs.end(); j++)
        {
          const Deriv& dTest2 = *j;
              
          MultipleDeriv m;
          m.put(dTest1);
          m.put(dTest2);
          if (e->hasNonzeroDeriv(m))
            {
              badDerivs.put(m);
            }
        }
    }
  
  TEST_FOR_EXCEPTION(badDerivs.size() != 0, RuntimeError,
                     "The expression " << expr << " is nonlinear in one "
                     "or more test functions. Nonlinear terms involving "
                     "the test functions " << endl << badDerivs << endl 
                     << " have been detected. Test functions must appear "
                     "only to the first degree");

  //  cerr << "max diff order=" << maxDiffOrder << endl;
  // cerr << "found derivs " << endl << nonzeroDerivs << endl;
  return nonzeroDerivs;
  
  
}
