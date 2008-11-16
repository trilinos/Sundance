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




DerivSet SymbPreprocessor::setupFwdProblem(const Expr& expr, 
                                           const Expr& tests,
                                           const Expr& unks,
                                           const Expr& unkEvalPts, 
                                           const Expr& unkParams,
                                           const Expr& unkParamEvalPts,
                                           const Expr& fixedParams,
                                           const Expr& fixedParamEvalPts,
                                           const Expr& fixedFields,
                                           const Expr& fixedFieldEvalPts,
                                           const EvalContext& context)
{
  Expr zero;
  Expr v = tests.flatten();
  Array<Expr> z(v.size());
  for (unsigned int i=0; i<v.size(); i++) z[i] = new ZeroExpr();
  zero = new ListExpr(z);

  return setupVariations(expr, 
                         tests, zero,
                         unks, unkEvalPts,
                         unkParams, unkParamEvalPts,
                         fixedFields, fixedFieldEvalPts,
                         fixedParams, fixedParamEvalPts,
                         context);
}




DerivSet SymbPreprocessor::setupSensitivities(const Expr& expr, 
                                              const Expr& tests,
                                              const Expr& unks,
                                              const Expr& unkEvalPts, 
                                              const Expr& unkParams,
                                              const Expr& unkParamEvalPts,
                                              const Expr& fixedParams,
                                              const Expr& fixedParamEvalPts,
                                              const Expr& fixedFields,
                                              const Expr& fixedFieldEvalPts,
                                              const EvalContext& context)
{
  Expr zero;
  Expr v = tests.flatten();
  Array<Expr> z(v.size());
  for (unsigned int i=0; i<v.size(); i++) z[i] = new ZeroExpr();
  zero = new ListExpr(z);

  return setupVariations(expr, 
                         tests, zero,
                         unks, unkEvalPts,
                         unkParams, unkParamEvalPts,
                         fixedFields, fixedFieldEvalPts,
                         fixedParams, fixedParamEvalPts,
                         context);
}


DerivSet SymbPreprocessor::setupFunctional(const Expr& expr, 
                                           const Expr& fixedParams,
                                           const Expr& fixedParamEvalPts,
                                           const Expr& fixedFields,
                                           const Expr& fixedFieldEvalPts,
                                           const EvalContext& context)
{
  Expr vars;
  Expr varEvalPts;
  Expr unks;
  Expr unkEvalPts;
  Expr unkParams;
  Expr unkParamEvalPts;

  return setupVariations(expr, 
                         vars, varEvalPts,
                         unks, unkEvalPts,
                         unkParams, unkParamEvalPts,
                         fixedFields, fixedFieldEvalPts,
                         fixedParams, fixedParamEvalPts,
                         context);
}





DerivSet SymbPreprocessor::setupGradient(const Expr& expr, 
                                         const Expr& vars,
                                         const Expr& varEvalPts,
                                         const Expr& fixedParams,
                                         const Expr& fixedParamEvalPts,
                                         const Expr& fixedFields,
                                         const Expr& fixedFieldEvalPts, 
                                         const EvalContext& context)
{
  Expr unks;
  Expr unkEvalPts;
  Expr unkParams;
  Expr unkParamEvalPts;

  return setupVariations(expr, 
                         vars, varEvalPts,
                         unks, unkEvalPts,
                         unkParams, unkParamEvalPts,
                         fixedFields, fixedFieldEvalPts,
                         fixedParams, fixedParamEvalPts,
                         context);
}




DerivSet SymbPreprocessor::setupVariations(const Expr& expr, 
                                           const Expr& vars,
                                           const Expr& varEvalPts,
                                           const Expr& unks,
                                           const Expr& unkEvalPts,
                                           const Expr& unkParams,
                                           const Expr& unkParamEvalPts,
                                           const Expr& fixedFields,
                                           const Expr& fixedFieldEvalPts, 
                                           const Expr& fixedParams,
                                           const Expr& fixedParamEvalPts, 
                                           const EvalContext& context)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  Array<Set<MultiSet<int> > > funcDerivs(3);
  Array<Set<MultiIndex> > spatialDerivs(3);

  SUNDANCE_OUT(Evaluator::classVerbosity() > VerbSilent,
               tab << "************ setting up variations of expr: " 
               << expr 
               << endl << tab << "context is " << context 
               << endl << tab << "vars are " << vars
               << endl << tab << "unks are " << unks
               << endl << tab << "unk parameters " << unkParams
               << endl << tab << "the eval points for the vars are " 
               << varEvalPts
               << endl << tab << "the eval points for the unks are " 
               << unkEvalPts
               << endl << tab 
               << " and the eval points for the parameters are " 
               << unkParamEvalPts);

  TEST_FOR_EXCEPTION(e==0, InternalError,
                     "Non-evaluatable expr " << expr.toString()
                     << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of variations, unknowns, parameters, and fixed fields */
  Expr v = vars.flatten();
  Expr v0 = varEvalPts.flatten();
  Expr u = unks.flatten();
  Expr u0 = unkEvalPts.flatten();
  Expr alpha = unkParams.flatten();
  Expr alpha0 = unkParamEvalPts.flatten();
  Expr beta = fixedParams.flatten();
  Expr beta0 = fixedParamEvalPts.flatten();
  Expr f = fixedFields.flatten();
  Expr f0 = fixedFieldEvalPts.flatten();
  

  SundanceUtils::Set<int> varID;
  SundanceUtils::Set<int> unkID;
  SundanceUtils::Set<int> unkParamID;
  SundanceUtils::Set<int> fixedParamID;
  SundanceUtils::Set<int> fixedID;


  /* check the var functions for redundancies and non-var funcs */
  for (unsigned int i=0; i<v.size(); i++)
    {
      const SymbolicFuncElement* vPtr
        = dynamic_cast<const SymbolicFuncElement*>(v[i].ptr().get());
      TEST_FOR_EXCEPTION(vPtr==0, RuntimeError, "list of variational funcs "
                         "contains a non-symbolic function " << v[i].toString());
      int fid = vPtr->funcComponentID();
      TEST_FOR_EXCEPTION(varID.contains(fid), RuntimeError,
                         "duplicate variational function in list "
                         << vars.toString());
      varID.put(fid);
      RefCountPtr<DiscreteFuncElement> v0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(v0[i].ptr());
      RefCountPtr<ZeroExpr> v0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(v0[i].ptr());
      TEST_FOR_EXCEPTION(v0Ptr.get()==NULL && v0ZeroPtr.get()==NULL,
                         RuntimeError,
                         "variation evaluation point " << u0[i].toString()
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
      int fid = uPtr->funcComponentID();
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
      int fid = aPtr->funcComponentID();
      TEST_FOR_EXCEPTION(unkParamID.contains(fid), RuntimeError,
                         "duplicate unknown parameter in list "
                         << alpha.toString());
      unkParamID.put(fid);
      RefCountPtr<Parameter> a0Ptr
        = rcp_dynamic_cast<Parameter>(alpha0[i].ptr());
      TEST_FOR_EXCEPTION(a0Ptr.get()==NULL,
                         RuntimeError,
                         "parameter evaluation point " << alpha0[i].toString()
                         << " is not a parameter");
      aPtr->substituteFunction(a0Ptr);
    }

  /* check the fixed parameters for redundancies and non-parameter funcs */
  for (unsigned int i=0; i<beta.size(); i++)
    {
      const UnknownParameterElement* bPtr
        = dynamic_cast<const UnknownParameterElement*>(beta[i].ptr().get());
      TEST_FOR_EXCEPTION(bPtr==0, RuntimeError,
                         "list of purported fixed unknown parameters "
                         "contains a function that is not an unknown parameter"
                         << beta[i].toString());
      int fid = bPtr->funcComponentID();
      TEST_FOR_EXCEPTION(fixedParamID.contains(fid), RuntimeError,
                         "duplicate fixed parameter in list "
                         << beta.toString());
      fixedParamID.put(fid);
      RefCountPtr<Parameter> b0Ptr
        = rcp_dynamic_cast<Parameter>(beta0[i].ptr());
      TEST_FOR_EXCEPTION(b0Ptr.get()==NULL,
                         RuntimeError,
                         "parameter evaluation point " << beta0[i].toString()
                         << " is not a parameter");
      bPtr->substituteFunction(b0Ptr);
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
      int fid = fPtr->funcComponentID();
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

  

  /* put together the set of functions that are active differentiation
   * variables */
  Array<SundanceUtils::Set<MultiSet<int> > > activeFuncIDs(3);
  activeFuncIDs[0].put(MultiSet<int>());
  if (context.topLevelDiffOrder() >= 1)
    {
      for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
        {
          activeFuncIDs[1].put(makeMultiSet<int>(*i));
          if (context.topLevelDiffOrder()==2)
            {
              for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
                {
                  activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
                }
              for (Set<int>::const_iterator 
                     j=unkParamID.begin(); j != unkParamID.end(); j++)
                {
                  activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
                }
            }
        }
    }
  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Finding nonzeros for expr " << endl);
  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Active funcs are " << activeFuncIDs << endl);

  Set<MultiIndex> miSet;
  miSet.put(MultiIndex());
  e->registerSpatialDerivs(context, miSet);
  
  Array<Set<MultipleDeriv> > RInput 
    = e->computeInputR(context, activeFuncIDs, spatialDerivs);

  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Top-level required funcs are " << RInput << endl);


  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Calling determineR()");
  
  e->determineR(context, RInput);


  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
               tab << endl << tab 
               << " ************* Setting up evaluators for expr " << endl);

  e->setupEval(context);

  if (verbosity<Evaluator>() > VerbMedium)
    { 
      SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
        tab << endl << tab 
        << " ************* Nonzeros are:");
      e->displayNonzeros(Out::os(), context);
    }

  DerivSet derivs = e->sparsitySuperset(context)->derivSet();


  SUNDANCE_OUT(verbosity<Evaluator>() > VerbLow,
    tab << endl << tab 
    << "Nonzero deriv set = " << derivs);

  return derivs;
}

