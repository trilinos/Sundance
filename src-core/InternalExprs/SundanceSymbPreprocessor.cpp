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

#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace Teuchos;





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
  const EvalContext& context,
  const ComputationType& compType)
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
    context, compType);
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
  const EvalContext& context,
  const ComputationType& compType)
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
    context, compType);
}


DerivSet SymbPreprocessor::setupFunctional(const Expr& expr, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts,
  const EvalContext& context,
  const ComputationType& compType)
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
    context,compType);
}





DerivSet SymbPreprocessor::setupGradient(const Expr& expr, 
  const Expr& vars,
  const Expr& varEvalPts,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts, 
  const EvalContext& context,
  const ComputationType& compType)
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
    context, compType);
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
  const EvalContext& context,
  const ComputationType& compType)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  Array<Set<MultiSet<int> > > funcDerivs(3);
  Array<Set<MultiIndex> > spatialDerivs(3);

  int verb=context.setupVerbosity();
  SUNDANCE_BANNER1(verb, tab, "in setupVariations()");
  verbosity<EvaluatableExpr>() = verb;
  SUNDANCE_MSG1(verb,
    tab << "************ setting up variations of expr: " 
    << expr 
    << endl << tab << "context is " << context 
    << endl << tab << "conp type is " << compType
    << endl << tab << "vars are " << vars
    << endl << tab << "unks are " << unks
    << endl << tab << "unk parameters " << unkParams
    << endl << tab << "fixed parameters " << fixedParams
    << endl << tab << "the eval points for the vars are " 
    << varEvalPts
    << endl << tab << "the eval points for the unks are " 
    << unkEvalPts
    << endl << tab 
    << "the eval points for the unknown parameters are " 
    << unkParamEvalPts 
    << endl << tab 
    << "the eval points for the fixed parameters are " 
    << fixedParamEvalPts 
    << tab << endl);

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
  


  Set<int> varID = processInputFuncs<SymbolicFuncElement>(v, v0);

  Set<int> unkID = processInputFuncs<UnknownFuncElement>(u, u0);

  Set<int> fixedID = processInputFuncs<UnknownFuncElement>(f, f0);

  Set<int> unkParamID 
    = processInputParams<UnknownParameterElement>(alpha, alpha0);

  Set<int> fixedParamID 
    = processInputParams<UnknownParameterElement>(beta, beta0);


  

  /* put together the set of functions that are active differentiation
   * variables */

  SUNDANCE_MSG5(verb, tab << "forming active set");
  Array<SundanceUtils::Set<MultiSet<int> > > activeFuncIDs(3);
  if (context.needsDerivOrder(0)) activeFuncIDs[0].put(MultiSet<int>());
  if (context.topLevelDiffOrder() >= 1)
  {
    for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      if (context.needsDerivOrder(1)) activeFuncIDs[1].put(makeMultiSet<int>(*i));
      if (context.topLevelDiffOrder()==2)
      {
        for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
        }
        if (compType==MatrixAndVector)
        {
          for (Set<int>::const_iterator 
                 j=unkParamID.begin(); j != unkParamID.end(); j++)
          {
            activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
          }
        }
        else if (compType==Sensitivities)
        {
          for (Set<int>::const_iterator 
                 j=fixedParamID.begin(); j != fixedParamID.end(); j++)
          {
            activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
          }
        }
      }
    }
  }
  SUNDANCE_MSG3(verb, tab << endl << tab 
    << " ************* Finding nonzeros for expr " << endl << tab);
  for (int i=0; i<=context.topLevelDiffOrder(); i++)
  {
    Tabs tab2;
    SUNDANCE_MSG4(verb, tab2 << "diff order=" << i << ", active funcs="
      << activeFuncIDs[i]);
  }

  Set<MultiIndex> miSet;
  miSet.put(MultiIndex());
  e->registerSpatialDerivs(context, miSet);
  
  SUNDANCE_MSG3(verb,
    tab << endl << tab 
    << " ************* finding required functions" << endl << tab);

  Array<Set<MultipleDeriv> > RInput 
    = e->computeInputR(context, activeFuncIDs, spatialDerivs);

  SUNDANCE_MSG3(verb,
    tab << endl << tab 
    << " ************* Top-level required funcs are " << RInput << endl << tab);


  SUNDANCE_MSG3(verb,
    tab << endl << tab 
    << " ************* Calling determineR()" << endl << tab);
  
  e->determineR(context, RInput);


  SUNDANCE_MSG3(verb,
    tab << endl << tab 
    << " ************* Setting up evaluators for expr " << endl << tab);

  e->setupEval(context);

  if (verb>1)
  { 
    SUNDANCE_MSG1(verb,
      tab << endl << tab 
      << " ************* Nonzeros are:");
    e->displayNonzeros(Out::os(), context);
  }

  DerivSet derivs = e->sparsitySuperset(context)->derivSet();


  SUNDANCE_MSG1(verb,
    tab << endl << tab 
    << "Nonzero deriv set = " << derivs);

  return derivs;
}


namespace SundanceCore
{
Expr makeZeros(const Expr& v)
{
  Array<Expr> z(v.size());
  for (unsigned int i=0; i<v.size(); i++) 
  {
    if (v[i].size()==1U) z[i] = new ZeroExpr();
    else z[i] = makeZeros(v[i]);
  }
  return new ListExpr(z);
}
}
