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
#include "SundanceUnknownParameterElement.hpp"
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
  Set<MultiSet<int> > activeFuncIDSet;
  activeFuncIDSet.put(MultiSet<int>());

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  e->findNonzeros(region, multiIndices, activeFuncIDSet, 
                  false);

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}

DerivSet SymbPreprocessor::setupFwdProblem(const Expr& expr, 
                                           const Expr& tests,
                                           const Expr& unks,
                                           const Expr& evalPts, 
                                           const Expr& unkParams,
                                           const Expr& unkParamEvalPts,
                                           const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               SUNDANCE_HEADER_LINE << tab 
               << "setting up forward problem expr: " << expr 
               << endl << tab << "with test functions " << tests
               << endl << tab << "and unknown functions " << unks
               << endl << tab 
               << "at the eval point " << evalPts
               << SUNDANCE_HEADER_LINE);


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
  Expr alpha = unkParams.flatten();
  Expr alpha0 = unkParamEvalPts.flatten();

  TEST_FOR_EXCEPTION(u.size() != u0.size(), RuntimeError,
                     "mismatched sizes of unknown list and evaluation points: "
                     "u.size() = " << u.size()
                     << ", u0.size() = " << u0.size());

  SundanceUtils::Set<int> testID;
  SundanceUtils::Set<int> unkID;
  SundanceUtils::Set<int> paramID;

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;

  /* check the test functions for redundancies and non-test funcs */
  for (unsigned int i=0; i<v.size(); i++)
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
  for (unsigned int i=0; i<u.size(); i++)
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


  /* check the unk parameters for redundancies and non-unk funcs */
  for (unsigned int i=0; i<alpha.size(); i++)
    {
      const UnknownParameterElement* aPtr
        = dynamic_cast<const UnknownParameterElement*>(alpha[i].ptr().get());
      TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
                         "list of purported unknown parameters "
                         "contains a function that is not an unknown parameter"
                         << alpha[i].toString());
      int fid = aPtr->funcID();
      TEST_FOR_EXCEPTION(paramID.contains(fid), RuntimeError,
                         "duplicate unknown parameter in list "
                         << alpha.toString());
      paramID.put(fid);
      RefCountPtr<Parameter> a0Ptr
        = rcp_dynamic_cast<Parameter>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL,
                         RuntimeError,
                         "parameter evaluation point " << alpha0[i].toString()
                         << " is not a parameter");
      aPtr->substituteFunction(a0Ptr);
    }

  for (Set<int>::const_iterator i=testID.begin(); i != testID.end(); i++)
    {
      MultiSet<int> testDeriv;
      testDeriv.put(*i);
      activeFuncIDs.put(testDeriv);
      for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          MultiSet<int> testUnkDeriv;
          testUnkDeriv.put(*i);
          testUnkDeriv.put(*j);
          activeFuncIDs.put(testUnkDeriv);
        }
    }

  

  /* Make sure there's no overlap between the test, unk, and params sets */
  for (Set<int>::const_iterator i=testID.begin(); i != testID.end(); i++)
    {
      TEST_FOR_EXCEPTION(unkID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and test lists");
      TEST_FOR_EXCEPTION(paramID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both test and param lists");
    }
  for (Set<int>::const_iterator i=unkID.begin(); i != unkID.end(); i++)
    {
      TEST_FOR_EXCEPTION(paramID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and param lists");
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

  SUNDANCE_OUT(verbosity<EvaluatableExpr>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  e->findNonzeros(region, multiIndices, activeFuncIDs, u0IsZero);

  
  SUNDANCE_OUT(verbosity<EvaluatableExpr>() > VerbLow,
               tab << endl << tab 
               << " ************* Sparsity pattern for expr " << endl
               << tab << e->toString() << endl
               << *(e->sparsitySuperset(region)) << endl
               << tab << " --------------------------------------------- "
               << endl);
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

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;
  activeFuncIDs.put(MultiSet<int>());


  /* check the unk functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<u.size(); i++)
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

  MultiSet<int> emptyDeriv;
  activeFuncIDs.put(emptyDeriv);

  for (Set<int>::const_iterator i=unkID.begin(); i != unkID.end(); i++)
    {
      MultiSet<int> unkDeriv;
      unkDeriv.put(*i);
      activeFuncIDs.put(unkDeriv);
      for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          MultiSet<int> unkUnkDeriv;
          unkUnkDeriv.put(*i);
          unkUnkDeriv.put(*j);
          activeFuncIDs.put(unkUnkDeriv);
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
  e->findNonzeros(region, multiIndices, activeFuncIDs, 
                  u0IsZero);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
  
}



DerivSet SymbPreprocessor::setupVariations(const Expr& expr, 
                                           const Expr& vars,
                                           const Expr& varEvalPts,
                                           const Expr& unks,
                                           const Expr& unkEvalPts, 
                                           const Expr& fixedFields,
                                           const Expr& fixedFieldEvalPts, 
                                           const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up variations of expr: " << expr 
               << endl << tab << "with respect to " << vars
               << endl << tab << ", solving for unknown functions " << unks
               << endl << tab << "the eval point for the unks is" 
               << unkEvalPts
               << endl << tab << " and the eval point for the variations is " 
               << varEvalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of variations, unknowns, and fixed fields */
  Expr v = vars.flatten();
  Expr u = unks.flatten();
  Expr alpha = fixedFields.flatten();
  Expr u0 = unkEvalPts.flatten();
  Expr v0 = varEvalPts.flatten();
  Expr alpha0 = fixedFieldEvalPts.flatten();

  SundanceUtils::Set<int> varID;
  SundanceUtils::Set<int> unkID;
  SundanceUtils::Set<int> fixedID;

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;

  /* check the variations for redundancies and non-unk funcs */
  for (unsigned int i=0; i<v.size(); i++)
    {
      const UnknownFuncElement* vPtr
        = dynamic_cast<const UnknownFuncElement*>(v[i].ptr().get());
      TEST_FOR_EXCEPTION(vPtr==0, RuntimeError, "list of variational funcs "
                         "contains a non-unknown function " 
                         << v[i].toString());
      int fid = vPtr->funcID();
      TEST_FOR_EXCEPTION(varID.contains(fid), RuntimeError,
                         "duplicate variation in list "
                         << vars.toString());
      varID.put(fid);
      RefCountPtr<DiscreteFuncElement> v0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(v0[i].ptr());
      RefCountPtr<ZeroExpr> v0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(v0[i].ptr());
      TEST_FOR_EXCEPTION(v0Ptr.get()==NULL && v0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "variational evaluation point " << v0[i].toString()
                         << " is neither a discrete function nor a zero expr");

      if (v0Ptr.get()==NULL)
        {
          vPtr->substituteZero();
        }
      else
        {
          vPtr->substituteFunction(v0Ptr);
        }
    }

  /* check the unk functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<u.size(); i++)
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
                         << u.toString());
      unkID.put(fid);
      RefCountPtr<DiscreteFuncElement> u0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(u0[i].ptr());
      RefCountPtr<ZeroExpr> u0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(u0[i].ptr());
      TEST_FOR_EXCEPTION(u0Ptr.get()==NULL && u0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "unknown evaluation point " << u0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (u0Ptr.get()==NULL)
        {
          uPtr->substituteZero();
        }
      else
        {
          uPtr->substituteFunction(u0Ptr);
        }
    }

  /* check the fixed functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<alpha.size(); i++)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(alpha[i].ptr().get());
      TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
                         "list of purported fixed funcs "
                         "contains a non-unknown function "
                         << alpha[i].toString());
      int fid = aPtr->funcID();
      TEST_FOR_EXCEPTION(fixedID.contains(fid), RuntimeError,
                         "duplicate unknown function in list "
                         << alpha.toString());
      fixedID.put(fid);
      RefCountPtr<DiscreteFuncElement> a0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(alpha0[i].ptr());
      RefCountPtr<ZeroExpr> a0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL && a0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "fixed-field evaluation point " 
                         << alpha0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (a0Ptr.get()==NULL)
        {
          aPtr->substituteZero();
        }
      else
        {
          aPtr->substituteFunction(a0Ptr);
        }
    }

  /* Make sure there's no overlap between the var, unk, and fixed sets */
  for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      TEST_FOR_EXCEPTION(unkID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and variation list");
    }
  for (Set<int>::const_iterator i=unkID.begin(); i != unkID.end(); i++)
    {
      TEST_FOR_EXCEPTION(varID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and variation list");
    }
  for (Set<int>::const_iterator i=unkID.begin(); i != unkID.end(); i++)
    {
      TEST_FOR_EXCEPTION(fixedID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and fixed list");
    }
  for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      TEST_FOR_EXCEPTION(fixedID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both variational and fixed list");
    }

  /* put together the set of functions that are active differentiation
   * variables */
  for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      MultiSet<int> varDeriv;
      varDeriv.put(*i);
      activeFuncIDs.put(varDeriv);
      for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          MultiSet<int> varUnkDeriv;
          varUnkDeriv.put(*i);
          varUnkDeriv.put(*j);
          activeFuncIDs.put(varUnkDeriv);
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
  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Active funcs are " << activeFuncIDs << endl);
  
  e->findNonzeros(region, multiIndices, activeFuncIDs, false);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}



DerivSet SymbPreprocessor::setupGradient(const Expr& expr, 
                                         const Expr& vars,
                                         const Expr& varEvalPts,
                                         const Expr& fixedFields,
                                         const Expr& fixedFieldEvalPts, 
                                         const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up gradient of expr: " << expr 
               << endl << tab << "with respect to " << vars
               << endl << tab << ". The value of the variable is " 
               << varEvalPts << " and the fields " << fixedFields 
               << " will be evaluated at " << fixedFieldEvalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of variations and fixed fields */
  Expr v = vars.flatten();
  Expr alpha = fixedFields.flatten();
  Expr v0 = varEvalPts.flatten();
  Expr alpha0 = fixedFieldEvalPts.flatten();

  SundanceUtils::Set<int> varID;
  SundanceUtils::Set<int> fixedID;

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;

  /* check the variations for redundancies and non-unk funcs */
  for (unsigned int i=0; i<v.size(); i++)
    {
      const UnknownFuncElement* vPtr
        = dynamic_cast<const UnknownFuncElement*>(v[i].ptr().get());
      TEST_FOR_EXCEPTION(vPtr==0, RuntimeError, "list of variational funcs "
                         "contains a non-unknown function " 
                         << v[i].toString());
      int fid = vPtr->funcID();
      TEST_FOR_EXCEPTION(varID.contains(fid), RuntimeError,
                         "duplicate variation in list "
                         << vars.toString());
      varID.put(fid);
      RefCountPtr<DiscreteFuncElement> v0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(v0[i].ptr());
      RefCountPtr<ZeroExpr> v0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(v0[i].ptr());
      TEST_FOR_EXCEPTION(v0Ptr.get()==NULL && v0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "variational evaluation point " << v0[i].toString()
                         << " is neither a discrete function nor a zero expr");

      if (v0Ptr.get()==NULL)
        {
          vPtr->substituteZero();
        }
      else
        {
          vPtr->substituteFunction(v0Ptr);
        }
    }

  /* check the fixed functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<alpha.size(); i++)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(alpha[i].ptr().get());
      TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
                         "list of purported fixed funcs "
                         "contains a non-unknown function "
                         << alpha[i].toString());
      int fid = aPtr->funcID();
      TEST_FOR_EXCEPTION(fixedID.contains(fid), RuntimeError,
                         "duplicate unknown function in list "
                         << alpha.toString());
      fixedID.put(fid);
      RefCountPtr<DiscreteFuncElement> a0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(alpha0[i].ptr());
      RefCountPtr<ZeroExpr> a0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL && a0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "fixed-field evaluation point " 
                         << alpha0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (a0Ptr.get()==NULL)
        {
          aPtr->substituteZero();
        }
      else
        {
          aPtr->substituteFunction(a0Ptr);
        }
    }

  /* Make sure there's no overlap between the var and fixed sets */
  for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      TEST_FOR_EXCEPTION(fixedID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both fixed and variation list");
    }

  /* put together the set of functions that are active differentiation
   * variables */
  MultiSet<int> zeroOrderDeriv;
  activeFuncIDs.put(zeroOrderDeriv);
  for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      MultiSet<int> varDeriv;
      varDeriv.put(*i);
      activeFuncIDs.put(varDeriv);
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
  e->findNonzeros(region, multiIndices, activeFuncIDs, false);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}



DerivSet SymbPreprocessor::setupFunctional(const Expr& expr, 
                                           const Expr& fixedFields,
                                           const Expr& fixedFieldEvalPts, 
                                           const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up functional: " << expr);

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent 
               && fixedFields.size()>0,
               "The fields " << fixedFields 
               << " will be evaluated at " << fixedFieldEvalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of fixed fields */
  Expr alpha = fixedFields.flatten();
  Expr alpha0 = fixedFieldEvalPts.flatten();

  SundanceUtils::Set<int> fixedID;

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;

  /* check the fixed functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<alpha.size(); i++)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(alpha[i].ptr().get());
      TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
                         "list of purported fixed funcs "
                         "contains a non-unknown function "
                         << alpha[i].toString());
      int fid = aPtr->funcID();
      TEST_FOR_EXCEPTION(fixedID.contains(fid), RuntimeError,
                         "duplicate unknown function in list "
                         << alpha.toString());
      fixedID.put(fid);
      RefCountPtr<DiscreteFuncElement> a0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(alpha0[i].ptr());
      RefCountPtr<ZeroExpr> a0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL && a0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "fixed-field evaluation point " 
                         << alpha0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (a0Ptr.get()==NULL)
        {
          aPtr->substituteZero();
        }
      else
        {
          aPtr->substituteFunction(a0Ptr);
        }
    }


  /* put together the set of functions that are active differentiation
   * variables */
  MultiSet<int> zeroOrderDeriv;
  activeFuncIDs.put(zeroOrderDeriv);
  

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
  e->findNonzeros(region, multiIndices, activeFuncIDs, false);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}







DerivSet SymbPreprocessor::setupSensitivities(const Expr& expr, 
                                              const Expr& tests,
                                              const Expr& unks,
                                              const Expr& unkEvalPts,
                                              const Expr& unkParams,
                                              const Expr& unkParamEvalPts,
                                              const Expr& fixedFields,
                                              const Expr& fixedFieldEvalPts, 
                                              const EvalContext& region)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up sensitivity calculation with expr: " 
               << expr 
               << endl << tab << "computing sensitivities of " << unks
               << endl << tab << "with respect to parameters " << unkParams
               << endl << tab << "the eval point for the unks is" 
               << unkEvalPts
               << endl << tab << " and the eval point for the parameters is " 
               << unkParamEvalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of variations, unknowns, and fixed fields */
  Expr v = tests.flatten();
  Expr u = unks.flatten();
  Expr u0 = unkEvalPts.flatten();
  Expr alpha = unkParams.flatten();
  Expr alpha0 = unkParamEvalPts.flatten();
  Expr f = fixedFields.flatten();
  Expr f0 = fixedFieldEvalPts.flatten();
  

  SundanceUtils::Set<int> testID;
  SundanceUtils::Set<int> unkID;
  SundanceUtils::Set<int> paramID;
  SundanceUtils::Set<int> fixedID;

  SundanceUtils::Set<MultiSet<int> > activeFuncIDs;

  /* check the test functions for redundancies and non-test funcs */
  for (unsigned int i=0; i<v.size(); i++)
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
  for (unsigned int i=0; i<u.size(); i++)
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
                         << u.toString());
      unkID.put(fid);
      RefCountPtr<DiscreteFuncElement> u0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(u0[i].ptr());
      RefCountPtr<ZeroExpr> u0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(u0[i].ptr());
      TEST_FOR_EXCEPTION(u0Ptr.get()==NULL && u0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "unknown evaluation point " << u0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (u0Ptr.get()==NULL)
        {
          uPtr->substituteZero();
        }
      else
        {
          uPtr->substituteFunction(u0Ptr);
        }
    }

  /* check the unk parameters for redundancies and non-unk funcs */
  for (unsigned int i=0; i<alpha.size(); i++)
    {
      const UnknownParameterElement* aPtr
        = dynamic_cast<const UnknownParameterElement*>(alpha[i].ptr().get());
      TEST_FOR_EXCEPTION(aPtr==0, RuntimeError,
                         "list of purported unknown parameters "
                         "contains a function that is not an unknown parameter"
                         << alpha[i].toString());
      int fid = aPtr->funcID();
      TEST_FOR_EXCEPTION(paramID.contains(fid), RuntimeError,
                         "duplicate unknown parameter in list "
                         << alpha.toString());
      paramID.put(fid);
      RefCountPtr<Parameter> a0Ptr
        = rcp_dynamic_cast<Parameter>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL,
                         RuntimeError,
                         "parameter evaluation point " << alpha0[i].toString()
                         << " is not a parameter");
      aPtr->substituteFunction(a0Ptr);
    }

  /* check the fixed functions for redundancies and non-unk funcs */
  for (unsigned int i=0; i<f.size(); i++)
    {
      const UnknownFuncElement* fPtr
        = dynamic_cast<const UnknownFuncElement*>(f[i].ptr().get());
      TEST_FOR_EXCEPTION(fPtr==0, RuntimeError,
                         "list of purported fixed funcs "
                         "contains a non-unknown function "
                         << f[i].toString());
      int fid = fPtr->funcID();
      TEST_FOR_EXCEPTION(fixedID.contains(fid), RuntimeError,
                         "duplicate unknown function in list "
                         << f.toString());
      fixedID.put(fid);
      RefCountPtr<DiscreteFuncElement> f0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(f0[i].ptr());
      RefCountPtr<ZeroExpr> f0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(f0[i].ptr());
      TEST_FOR_EXCEPTION(f0Ptr.get()==NULL && f0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "fixed-field evaluation point " 
                         << f0[i].toString()
                         << " is neither a discrete function nor a zero expr");
      if (f0Ptr.get()==NULL)
        {
          fPtr->substituteZero();
        }
      else
        {
          fPtr->substituteFunction(f0Ptr);
        }
    }

  /* Make sure there's no overlap between the test, unk, params, and fixed sets */
  for (Set<int>::const_iterator i=testID.begin(); i != testID.end(); i++)
    {
      TEST_FOR_EXCEPTION(unkID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and test lists");
      TEST_FOR_EXCEPTION(paramID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both test and param lists");
      TEST_FOR_EXCEPTION(fixedID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both test and fixed lists");
    }
  for (Set<int>::const_iterator i=unkID.begin(); i != unkID.end(); i++)
    {
      TEST_FOR_EXCEPTION(fixedID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and fixed lists");
      TEST_FOR_EXCEPTION(paramID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both unknown and param lists");
    }
  for (Set<int>::const_iterator i=fixedID.begin(); i != fixedID.end(); i++)
    {
      TEST_FOR_EXCEPTION(paramID.contains(*i), RuntimeError,
                         "Function with ID=" << *i << " appears in "
                         "both param and fixed lists");
    }
  

  /* put together the set of functions that are active differentiation
   * variables */
  for (Set<int>::const_iterator i=testID.begin(); i != testID.end(); i++)
    {
      MultiSet<int> testDeriv;
      for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          MultiSet<int> testUnkDeriv;
          testUnkDeriv.put(*i);
          testUnkDeriv.put(*j);
          activeFuncIDs.put(testUnkDeriv);
        }
      for (Set<int>::const_iterator j=paramID.begin(); j != paramID.end(); j++)
        {
          MultiSet<int> testUnkDeriv;
          testUnkDeriv.put(*i);
          testUnkDeriv.put(*j);
          activeFuncIDs.put(testUnkDeriv);
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
  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Active funcs are " << activeFuncIDs << endl);
  
  e->findNonzeros(region, multiIndices, activeFuncIDs, false);

   SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
                << " ************* Setting up evaluators for expr " << endl);
  e->setupEval(region);

  DerivSet derivs = e->sparsitySuperset(region)->derivSet();

  return derivs;
}

