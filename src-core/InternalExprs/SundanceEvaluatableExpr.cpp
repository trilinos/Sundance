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

#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSymbolicFuncElement.hpp"
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




EvaluatableExpr::EvaluatableExpr()
	: ScalarExpr(), 
    evaluators_(),
    sparsity_(),
    orderOfDependency_(MultiIndex::maxDim(), -1),
    knownNonzeros_(),
    funcIDSet_(),
    nodesHaveBeenCounted_(false)
{
  addFuncIDCombo(MultiSet<int>());
}




RefCountPtr<SparsitySubset> 
EvaluatableExpr::sparsitySubset(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs) const 
{
  Tabs tab;
  RefCountPtr<SparsitySuperset> super = sparsitySuperset(context);
  SUNDANCE_VERB_HIGH(tab << "getting subset for miSet=" 
                     << multiIndices.toString() << " and active funcs="
                     << activeFuncIDs << ". Superset is "
                     << *super);
  if (!super->hasSubset(multiIndices, activeFuncIDs))
    {
      SUNDANCE_VERB_HIGH(tab << "miSet, func combination not found");
      super->addSubset(multiIndices, activeFuncIDs);
    }


  return super->subset(multiIndices, activeFuncIDs);
}




RefCountPtr<SparsitySuperset> 
EvaluatableExpr::sparsitySuperset(const EvalContext& context) const 
{
  RefCountPtr<SparsitySuperset> rtn;

  if (sparsity_.containsKey(context))
    {
      rtn = sparsity_.get(context);
    }
  else
    {
      rtn = rcp(new SparsitySuperset());
      sparsity_.put(context, rtn);
    }
  return rtn;
}



 const RefCountPtr<Evaluator>&
EvaluatableExpr::evaluator(const EvalContext& context) const 
{
  TEST_FOR_EXCEPTION(!evaluators_.containsKey(context), RuntimeError, 
                     "Evaluator not found for context " << context);
  return evaluators_.get(context);
}

bool EvaluatableExpr
::nonzerosAreKnown(const EvalContext& context,
                   const Set<MultiIndex>& multiIndices,
                   const Set<MultiSet<int> >& activeFuncIDs,
                   bool regardFuncsAsConstant) const 
{
  NonzeroSpecifier spec(context, multiIndices, 
                        activeFuncIDs, regardFuncsAsConstant);
  return knownNonzeros_.contains(spec);
}

void EvaluatableExpr::addKnownNonzero(const EvalContext& context,
                                      const Set<MultiIndex>& multiIndices,
                                      const Set<MultiSet<int> >& activeFuncIDs,
                                      bool regardFuncsAsConstant) const 
{
  NonzeroSpecifier spec(context, multiIndices, 
                        activeFuncIDs, regardFuncsAsConstant);
  knownNonzeros_.put(spec);
}

void EvaluatableExpr::evaluate(const EvalManager& mgr,
                               Array<double>& constantResults,
                               Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  evaluator(mgr.getRegion())->eval(mgr, constantResults, vectorResults);
}

void EvaluatableExpr::setupEval(const EvalContext& context) const
{
  if (!evaluators_.containsKey(context))
    {
      RefCountPtr<Evaluator> eval = rcp(createEvaluator(this, context));
      evaluators_.put(context, eval);
    }
}

void EvaluatableExpr::showSparsity(ostream& os, 
                                   const EvalContext& context) const
{
  Tabs tab0;
  os << tab0 << "Node: " << toString() << endl;
  sparsitySuperset(context)->displayAll(os);
}



int EvaluatableExpr::maxOrder(const Set<MultiIndex>& m) const 
{
  int rtn = 0;
  for (Set<MultiIndex>::const_iterator i=m.begin(); i != m.end(); i++)
    {
      rtn = max(i->order(), rtn);
    }
  return rtn;
}

void EvaluatableExpr::addFuncIDCombo(const MultiSet<int>& funcIDSet)
{
  funcIDSet_.put(funcIDSet);
  for (MultiSet<int>::const_iterator 
         i=funcIDSet.begin(); i != funcIDSet.end(); i++)
    {
      funcDependencies_.put(*i);
    }
}


void EvaluatableExpr::setFuncIDSet(const Set<MultiSet<int> >& funcIDSet)
{
  funcIDSet_ = funcIDSet;
  for (Set<MultiSet<int> >::const_iterator 
         i=funcIDSet.begin(); i != funcIDSet.begin(); i++)
    {
      const MultiSet<int>& d = *i;
      for (MultiSet<int>::const_iterator j=d.begin(); j != d.begin(); j++)
        {
          funcDependencies_.put(*j);
        }
    }
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


bool EvaluatableExpr::isEvaluatable(const ExprBase* expr) 
{
  return dynamic_cast<const EvaluatableExpr*>(expr) != 0;
}


int EvaluatableExpr::countNodes() const
{
  nodesHaveBeenCounted_ = true;
  return 1;
}


RefCountPtr<Set<int> > EvaluatableExpr::getFuncIDSet(const Expr& funcs)
{
  RefCountPtr<Set<int> > rtn = rcp(new Set<int>());

  Expr f = funcs.flatten();
  for (int i=0; i<f.size(); i++)
    {
      const SymbolicFuncElement* sfe 
        = dynamic_cast<const SymbolicFuncElement*>(f[i].ptr().get());
      TEST_FOR_EXCEPTION(sfe==0, RuntimeError,
                         "non-symbolic function expr " << f[i]
                         << " found in getFuncIDSet()");
      rtn->put(sfe->funcID());
    }
  return rtn;
}



