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

#include "SundanceExprWithChildren.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceNullEvaluator.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnaryExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;



ExprWithChildren::ExprWithChildren(const Array<RefCountPtr<ScalarExpr> >& children)
	: EvaluatableExpr(), 
    children_(children),
    contextToQWMap_(4)
{}


bool ExprWithChildren::isConstant() const
{
  for (unsigned int i=0; i<children_.size(); i++) 
    {
      if (!children_[i]->isConstant()) return false;
    }
  return true;
}

void ExprWithChildren::accumulateFuncSet(Set<int>& funcIDs,
                                         const Set<int>& activeFuncs) const
{
  for (unsigned int i=0; i<children_.size(); i++) 
    {
      children_[i]->accumulateFuncSet(funcIDs, activeFuncs);
    }
}


const EvaluatableExpr* ExprWithChildren::evaluatableChild(int i) const
{
  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(children_[i].get());

  TEST_FOR_EXCEPTION(e==0, InternalError, 
                     "ExprWithChildren: cast of child [" 
                     << children_[i]->toString()
                     << " to evaluatable expr failed");

  return e;
}


void ExprWithChildren::findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    const Set<MultiSet<int> >& inputActiveFuncIDs,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;
  ExprWithChildren* nonConstSelf = const_cast<ExprWithChildren*>(this);
  EvaluatableExpr* nonConstBase = static_cast<EvaluatableExpr*>(nonConstSelf);
  nonConstBase->verbosity() = TSFExtended::verbosity<EvaluatableExpr>();
  SUNDANCE_VERB_HIGH(tabs << "finding nonzeros for ExprWithChildren " 
                       << toString() << " subject to multiindex set "
                       << multiIndices.toString());
  if (this->verbosity() > VerbHigh)
    {
      SUNDANCE_OUT(true, tabs << "num children = " << children_.size());
      for (unsigned int i=0; i<children_.size(); i++)
        {
          Tabs tab1;
          SUNDANCE_OUT(true, tab1 << "child #" << i << " = " 
                       << evaluatableChild(i)->toString());
        }
    }


  Set<MultiSet<int> > filteredActiveFuncs = filterActiveFuncs(inputActiveFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, filteredActiveFuncs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_HIGH(tabs << "...reusing previously computed data");
      return;
    }

  const UnaryExpr* ue = dynamic_cast<const UnaryExpr*>(this);
  if (ue != 0) ue->addActiveFuncs(context, filteredActiveFuncs);

  RefCountPtr<SparsitySubset> subset 
    = sparsitySubset(context, multiIndices, filteredActiveFuncs, false);

  /* The sparsity pattern is the union of the 
   * operands' sparsity patterns. If any functional derivatives
   * appear in multiple operands, the state of that derivative is
   * the more general of the states */
  for (unsigned int i=0; i<children_.size(); i++)
    {
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "finding nonzeros for child " 
                           << evaluatableChild(i)->toString());
      evaluatableChild(i)->findNonzeros(context, multiIndices,
                                        filteredActiveFuncs,
                                        regardFuncsAsConstant);

      Set<MultiSet<int> > filteredChildActiveFuncs 
        = evaluatableChild(i)->filterActiveFuncs(filteredActiveFuncs);
      
      SUNDANCE_VERB_HIGH(tab1 << "reading sparsity subset for child " 
                           << evaluatableChild(i)->toString());

      RefCountPtr<SparsitySubset> childSparsitySubset 
        = evaluatableChild(i)->sparsitySubset(context, multiIndices, filteredChildActiveFuncs, true);
          

      SUNDANCE_VERB_HIGH(tabs << "child #" << i 
                           << " sparsity subset is " 
                           << endl << *childSparsitySubset);
      
      for (int j=0; j<childSparsitySubset->numDerivs(); j++)
        {
          subset->addDeriv(childSparsitySubset->deriv(j),
                           childSparsitySubset->state(j));
        }
    }

  SUNDANCE_VERB_HIGH(tabs << "done with ExprWithChildren " 
                     + toString() << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "done with ExprWithChildren " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, filteredActiveFuncs,
                  regardFuncsAsConstant);
}

void ExprWithChildren::setupEval(const EvalContext& context) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "expr " + toString() 
                     << ": creating evaluators for children");
  SUNDANCE_VERB_HIGH(tabs << "my sparsity superset = " << *sparsitySuperset(context));


  if (sparsitySuperset(context)->numDerivs() > 0)
    {
      for (unsigned int i=0; i<children_.size(); i++)
        {
          Tabs tabs1;
          SUNDANCE_VERB_HIGH(tabs1 << "creating evaluator for child " 
                             << evaluatableChild(i)->toString());
          evaluatableChild(i)->setupEval(context);
        }
    }

  if (!evaluators().containsKey(context))
    {
      RefCountPtr<Evaluator> eval;
      if (sparsitySuperset(context)->numDerivs()>0)
        {
          eval = rcp(createEvaluator(this, context));
        }
      else
        {
          eval = rcp(new NullEvaluator());
        }
      evaluators().put(context, eval);
    }
}

void ExprWithChildren::showSparsity(ostream& os, 
                                    const EvalContext& context) const
{
  Tabs tab0;
  os << tab0 << "Node: " << toString() << endl;
  sparsitySuperset(context)->displayAll(os);
  for (unsigned int i=0; i<children_.size(); i++)
    {
      Tabs tab1;
      os << tab1 << "Child " << i << endl;
      evaluatableChild(i)->showSparsity(os, context);
    }
}


bool ExprWithChildren::allTermsHaveTestFunctions() const
{
  for (unsigned int i=0; i<children_.size(); i++)
    {
      if (evaluatableChild(i)->allTermsHaveTestFunctions()) return true;
    }
  return false;
}


bool ExprWithChildren::hasTestFunctions() const
{
  for (unsigned int i=0; i<children_.size(); i++)
    {
      if (evaluatableChild(i)->hasTestFunctions()) return true;
    }
  return false;
}

void ExprWithChildren::getUnknowns(Set<int>& unkID, Array<Expr>& unks) const
{
  for (unsigned int i=0; i<children_.size(); i++)
    {
      const RefCountPtr<ExprBase>& e = children_[i];
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(e.get());
      if (u != 0)
        {
          Expr expr(e);
          if (!unkID.contains(u->funcID())) 
            {
              unks.append(expr);
              unkID.put(u->funcID());
            }
        }
      evaluatableChild(i)->getUnknowns(unkID, unks);
    }
  
}


int ExprWithChildren::countNodes() const
{
  if (nodesHaveBeenCounted()) 
    {
      return 0;
    }

  /* count self */
  int count = EvaluatableExpr::countNodes();

  /* count children */
  for (unsigned int i=0; i<children_.size(); i++)
    {
      if (!evaluatableChild(i)->nodesHaveBeenCounted())
        {
          count += evaluatableChild(i)->countNodes();
        }
    }
  return count;
}


Set<MultipleDeriv> 
ExprWithChildren::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  
  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
    {
      /* If I am an arbitrary nonlinear expression, I cannot be known to
       * be zero regardless of the state of my arguments. Return the
       * zeroth-order derivative. */
      if (!(isLinear() || isProduct())) 
        {
          rtn.put(MultipleDeriv());
          return rtn;
        }

      /* At this point, I've dealt with arbitrary nonlinear exprs so 
       * I know I'm either a product or a linear expr */

      const Set<Array<int> >& Q = findQ_W(0, context);

      /* If there are no nonzero terms, a linear combination or a product
       * will be zero. Return the empty set. */
      if (Q.size()==0)
        {
          return rtn;
        }
      
      /* if I'm a linear combination and any term is nonzero, I am nonzero */
      if (isLinear())
        {
          rtn.put(MultipleDeriv());          
          return rtn;
        }

      /* The only possibility remaining is that I'm a product.
       * If any term is zero, I am zero. If a term is 
       * known to be zero, it's index will not appear in Q, so that 
       * comparing the size of Q and the number of children tells me
       * if I'm zero.
       */
      if (((int) Q.size()) == numChildren())
        {
          rtn.put(MultipleDeriv());
        }
      return rtn;
    }


  /* now do arbitrary order deriv */
  for (int i=1; i<=order; i++) 
    {
      const Set<Array<int> >& Q = findQ_W(i, context);
      
      for (Set<Array<int> >::const_iterator j=Q.begin(); j!=Q.end(); j++)
        {
          if (j->size()==1)
            {
              int childIndex = (*j)[0];
              const Set<MultipleDeriv>& childW 
                = evaluatableChild(childIndex)->findW(order, context);
              rtn.merge(childW);
            }
          else if (j->size()==2)
            {
              int i1 = (*j)[0];
              int i2 = (*j)[1];
              
              if (order==2)
                {
                  const Set<MultipleDeriv>& childW1 
                    = evaluatableChild(i1)->findW(1, context);
                  const Set<MultipleDeriv>& childW2 
                    = evaluatableChild(i2)->findW(1, context);
                  rtn.merge(setProduct(childW1, childW2));
                }
              else
                {
                  const Set<MultipleDeriv>& childW12 
                    = evaluatableChild(i1)->findW(2, context);
                  const Set<MultipleDeriv>& childW21 
                    = evaluatableChild(i2)->findW(1, context);
                  const Set<MultipleDeriv>& childW11 
                    = evaluatableChild(i1)->findW(1, context);
                  const Set<MultipleDeriv>& childW22 
                    = evaluatableChild(i2)->findW(2, context);
                  rtn.merge(setProduct(childW12, childW21));
                  rtn.merge(setProduct(childW11, childW22));
                }
            }
          else if (j->size()==3)
            {
              int i1 = (*j)[0];
              int i2 = (*j)[1];
              int i3 = (*j)[2];
              
              if (order==3)
                {
                  const Set<MultipleDeriv>& childW1 
                    = evaluatableChild(i1)->findW(1, context);
                  const Set<MultipleDeriv>& childW2 
                    = evaluatableChild(i2)->findW(1, context);
                  const Set<MultipleDeriv>& childW3 
                    = evaluatableChild(i3)->findW(1, context);
                  rtn.merge(setProduct(setProduct(childW1, childW2), childW3));
                }
              
            }
            
        }
    }

  return rtn;
}

const Set<Array<int> >& 
ExprWithChildren::findQ_W(int order, 
                          const EvalContext& context) const
{
  if (!contextToQWMap_[order].containsKey(context))
    {
      contextToQWMap_[order].put(context, internalFindQ_W(order, context));
    }
  return contextToQWMap_[order].get(context);
}

Set<Array<int> > ExprWithChildren
::internalFindQ_W(int order, 
                  const EvalContext& context) const
{
  Set<Array<int> > rtn;

  if (order <= 1)
    {
      for (int i=0; i<numChildren(); i++) rtn.put(tuple(i));
      return rtn;
    }
  if (order==2)
    {
      for (int i=0; i<numChildren(); i++) 
        {
          for (int j=0; j<=i; j++) 
            {
              rtn.put(tuple(i,j));
            }
        }
      return rtn;
    }
  if (order==3)
    {
      for (int i=0; i<numChildren(); i++) 
        {
          for (int j=0; j<=i; j++) 
            {
              for (int k=0; k<=j; k++) 
                {
                  rtn.put(tuple(i,j, k));
                }
            }
        }
      return rtn;
    }
  
  TEST_FOR_EXCEPTION(order > 3, RuntimeError, "invalid order " << order);
  
  return rtn;
}
