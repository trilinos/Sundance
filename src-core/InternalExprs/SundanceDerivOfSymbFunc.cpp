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

#include "SundanceDerivOfSymbFunc.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

DerivOfSymbFunc::DerivOfSymbFunc(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg)
  : DiffOp(op, arg), funcID_(-1)
{
  const SymbolicFuncElement* f 
    = dynamic_cast<const SymbolicFuncElement*>(evaluatableArg());
  TEST_FOR_EXCEPTION(f==0, InternalError, "argument to DerivOfSymbFunc ctor "
                     "is not a symbolic function");
  funcID_ = f->funcID();
}


FunctionalDeriv* DerivOfSymbFunc::representMeAsFunctionalDeriv() const 
{
  const FuncElementBase* f 
    = dynamic_cast<const FuncElementBase*>(evaluatableArg());
  TEST_FOR_EXCEPTION(f==0, InternalError, "DerivOfSymbFunc::"
                     "representMeAsFunctionalDeriv(), 'this' pointer "
                     "is not a symbolic function");
  return new FunctionalDeriv(f, mi());
}

Evaluator* DerivOfSymbFunc::createEvaluator(const EvaluatableExpr* expr,
                                   const EvalContext& context) const
{
  return new DerivOfSymbFuncEvaluator(dynamic_cast<const DerivOfSymbFunc*>(expr), context);
}

Set<MultiSet<int> > DerivOfSymbFunc
::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs) const
{
  Tabs tabs;
  Set<MultiSet<int> > rtn;

  SUNDANCE_VERB_MEDIUM(tabs << "arg dependencies are " 
                       << evaluatableArg()->funcDependencies().toString());

  for (Set<int>::const_iterator 
         j=evaluatableArg()->funcDependencies().begin(); 
       j != evaluatableArg()->funcDependencies().end(); j++)
    {
      MultiSet<int> f;
      f.put(*j);
      rtn.put(f);
    }
  
  for (Set<MultiSet<int> >::const_iterator 
         i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& d = *i;
      if (d.size() >= 1) continue;
      SUNDANCE_VERB_MEDIUM(tabs << "deriv from outside is " << d.toString());
      rtn.put(d);
    }
  return rtn;
}

void DerivOfSymbFunc::findNonzeros(const EvalContext& context,
                          const Set<MultiIndex>& multiIndices,
                          const Set<MultiSet<int> >& activeFuncIDs,
                          bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for deriv of func " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());
  SUNDANCE_VERB_MEDIUM(tabs << "active funcs are " << activeFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  addActiveFuncs(context, activeFuncIDs);
  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs);


  /* Figure out the sparsity pattern for the argument.  */
  Set<MultiIndex> argMI = argMultiIndices(multiIndices);
  
  SUNDANCE_VERB_MEDIUM(tabs << "arg multi index set is " << endl << argMI);

  

  Set<MultiSet<int> > argFuncs = argActiveFuncs(activeFuncIDs);

  SUNDANCE_VERB_MEDIUM(tabs << "arg active func set is " << endl << argFuncs);

  evaluatableArg()->findNonzeros(context, argMI,
                                 argFuncs,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsity
    = evaluatableArg()->sparsitySubset(context, argMI, argFuncs);

  SUNDANCE_VERB_MEDIUM(tabs << "arg sparsity subset is " 
                       << endl << *argSparsity);

  for (int i=0; i<argSparsity->numDerivs(); i++)
    {
      Tabs tab1;

      const MultipleDeriv& md = argSparsity->deriv(i);

      if (md.order()==0) continue;
      
      SUNDANCE_VERB_MEDIUM(tab1 << "finding the effect of the argument's "
                           "nonzero derivative " << md);
      TEST_FOR_EXCEPTION(md.order() > 1, InternalError,
                         "detected a functional derivative of order > 1 "
                         "acting on a symbolic function. This shouldn't "
                         "happen");
      const Deriv& d = *(md.begin());
      TEST_FOR_EXCEPTION(d.isCoordDeriv(), InternalError, 
                         "CoordDeriv detected in DerivOfSymbFunc");
      const FuncElementBase* f = d.funcDeriv()->func();
      const SymbolicFuncElement* s = dynamic_cast<const SymbolicFuncElement*>(f);
      TEST_FOR_EXCEPTION(s==0, InternalError, "argument of DerivOfSymbFunc "
                         "is not a symbolic func");
      Deriv dMi = d.funcDeriv()->derivWrtMultiIndex(mi());
      subset->addDeriv(dMi, ConstantDeriv);
      if (!s->evalPtIsZero())
        {
          subset->addDeriv(MultipleDeriv(), VectorDeriv);
        }
    }


  SUNDANCE_VERB_HIGH(tabs << "DerivOfSymbFunc " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "DerivOfSymbFunc " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}

