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
using namespace TSFExtended;



DerivSet SymbPreprocessor::setupExpr(const Expr& expr, 
                                     const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");
  Set<MultiIndex> multiIndices;
  multiIndices.put(MultiIndex());

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  e->findNonzeros(region, multiIndices, false);

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}

DerivSet SymbPreprocessor::setupExpr(const Expr& expr, 
                                     const Expr& tests,
                                     const Expr& unks,
                                     const Expr& evalPts, 
                                     const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up expr: " << expr 
               << endl << tab << "with test functions " << tests
               << endl << tab << "and unknown functions " << unks
               << endl << tab << "at the eval point " << evalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  bool missingTestFunction = !e->allTermsHaveTestFunctions();
  TEST_FOR_EXCEPTION(missingTestFunction, RuntimeError,
                     "Weak form " << expr << " contains at least one "
                     "term without a test function");

  bool u0IsZero = false;

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

  /* ----------------------------------------------------------
   * set up the expression for evaluation 
   * ----------------------------------------------------------
   *
   * In this step, we work out the sparsity structure of the expression's
   * functional derivatives, and create evaluator objects. 
   */


  
  Set<MultiIndex> multiIndices;
  multiIndices.put(MultiIndex());

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  e->findNonzeros(region, multiIndices, u0IsZero);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}


DerivSet SymbPreprocessor::setupExpr(const Expr& expr, 
                                     const Expr& unks,
                                     const Expr& evalPts, 
                                     const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  bool u0IsZero = false;

  /* make flat lists of unknowns */
  Expr u = unks.flatten();
  Expr u0 = evalPts.flatten();

  SundanceUtils::Set<int> unkID;

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

  /* ----------------------------------------------------------
   * set up the expression for evaluation 
   * ----------------------------------------------------------
   *
   * In this step, we work out the sparsity structure of the expression's
   * functional derivatives, and create evaluator objects. 
   */


  
  Set<MultiIndex> multiIndices;
  multiIndices.put(MultiIndex());

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  e->findNonzeros(region, multiIndices, u0IsZero);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
  
}



