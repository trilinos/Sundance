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
#include "SundanceCombinatorialUtils.hpp"
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
    contextToQWMap_(4),
    contextToQVMap_(4),
    contextToQCMap_(4)
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
ExprWithChildren::product(const Array<int>& J, const Array<int>& K,
                          DerivSubsetSpecifier dss,
                          const EvalContext& context) const
{
  TEST_FOR_EXCEPTION(J.size() != K.size(), InternalError,
                     "mismatched index set sizes");
  TEST_FOR_EXCEPTION(J.size() == 0, InternalError,
                     "empty index set");

  Set<MultipleDeriv> rtn 
    = evaluatableChild(J[0])->findDerivSubset(K[0], dss, context);
  
  for (unsigned int i=1; i<J.size(); i++)
    {
      const Set<MultipleDeriv>& S 
        = evaluatableChild(J[i])->findDerivSubset(K[i], dss, context);
      rtn = setProduct(rtn, S);
    }

  return rtn;
}

Set<MultipleDeriv> 
ExprWithChildren::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "EWC::internalFindV() for " << toString());

  Set<MultipleDeriv> rtn;



  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
    {
      for (int i=0; i<numChildren(); i++)
        {
          if (!childIsRequired(i, order, context)) continue;
          const Set<MultipleDeriv>& childV0 
            = evaluatableChild(i)->findV(0, context);        
          rtn.merge(childV0);
        }
      const Set<MultipleDeriv>& R0 = findR(0, context);
      rtn = rtn.intersection(R0);

      SUNDANCE_VERB_HIGH(tab0 << "V[" << order << "]=" << rtn);
      SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindV for "
                         << toString());
      return rtn;
    }


  /* now do arbitrary order derivatives with the multivariable chain rule*/
  Array<Array<Array<int> > > comp = compositions(order);
  for (int i=1; i<=order; i++) 
    {
      Tabs tab1;
      SUNDANCE_VERB_EXTREME(tab1 << "i=" << i);
      const Set<MultiSet<int> >& QC = findQ_C(i, context);
      SUNDANCE_VERB_EXTREME(tab1 << "QC = " << QC);
      for (Set<MultiSet<int> >::const_iterator j=QC.begin(); j!=QC.end(); j++)
        {
          Tabs tab2;
          Array<int> J = j->elements();
          const Array<Array<int> >& K = comp[J.size()-1];
          SUNDANCE_VERB_EXTREME(tab2 << "J=" << J);
          SUNDANCE_VERB_EXTREME(tab2 << "K=" << K);

          for (unsigned int k=0; k<K.size(); k++)
            {
              Tabs tab3;
              Set<MultipleDeriv> RProd = product(J, K[k], 
                                                 RequiredNonzeros, context);
              Set<MultipleDeriv> CProd = product(J, K[k], 
                                                 ConstantNonzeros, context);
              SUNDANCE_VERB_EXTREME(tab3 << "CProd = " << CProd);
              SUNDANCE_VERB_EXTREME(tab3 << "RProd = " << RProd);
              rtn.merge(RProd.setDifference(CProd));
            }
        }

      const Set<MultiSet<int> >& QV = findQ_V(i, context);
      SUNDANCE_VERB_EXTREME(tab1 << "QV = " << QV);
      for (Set<MultiSet<int> >::const_iterator j=QV.begin(); j!=QV.end(); j++)
        {
          Tabs tab2;
          Array<int> J = j->elements();
          const Array<Array<int> >& K = comp[J.size()-1];
          SUNDANCE_VERB_EXTREME(tab2 << "J=" << J);
          SUNDANCE_VERB_EXTREME(tab2 << "K=" << K);

          for (unsigned int k=0; k<K.size(); k++)
            {
              Tabs tab3;
              Set<MultipleDeriv> RProd = product(J, K[k], 
                                                 RequiredNonzeros, context);
              SUNDANCE_VERB_EXTREME(tab3 << "RProd = " << RProd);
              rtn.merge(RProd);
            }
        }
    }

  SUNDANCE_VERB_HIGH(tab0 << "V[" << order << "]=" << rtn);
  SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindV for "
                     << toString());

  return rtn;
}


Set<MultipleDeriv> 
ExprWithChildren::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab0;
  Set<MultipleDeriv> rtn;
  bool started = false;

  SUNDANCE_VERB_HIGH(tab0 << "EWC::internalFindC() for " << toString());

  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
    {
      for (int i=0; i<numChildren(); i++)
        {
          if (!childIsRequired(i, 0, context)) continue;
          const Set<MultipleDeriv>& childV0 
            = evaluatableChild(i)->findV(0, context);        
          rtn.merge(childV0);
        }
      const Set<MultipleDeriv>& R0 = findR(0, context);
      rtn = R0.setDifference(rtn);
      SUNDANCE_VERB_HIGH(tab0 << "C[" << order << "]=" << rtn);
      SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindC for "
                         << toString());
      return rtn;
    }


  /* now do arbitrary order derivatives with the multivariable chain rule*/
  Array<Array<Array<int> > > comp = compositions(order);
  const Set<MultipleDeriv>& RM = findR(order, context);
  for (int i=1; i<=order; i++) 
    {
      Tabs tab1;
      SUNDANCE_VERB_EXTREME(tab1 << "i=" << i);


      const Set<MultiSet<int> >& QC = findQ_C(i, context);
      SUNDANCE_VERB_EXTREME(tab1 << "finding CProd union (R\\RProd) over QC");      
      SUNDANCE_VERB_EXTREME(tab1 << "QC = " << QC);


      for (Set<MultiSet<int> >::const_iterator j=QC.begin(); j!=QC.end(); j++)
        {
          Tabs tab2;
          Array<int> J = j->elements();
          const Array<Array<int> >& K = comp[J.size()-1];
          SUNDANCE_VERB_EXTREME(tab2 << "J=" << J);
          SUNDANCE_VERB_EXTREME(tab2 << "K=" << K);

          for (unsigned int k=0; k<K.size(); k++)
            {
              Tabs tab3;
              Set<MultipleDeriv> RProd = product(J, K[k], 
                                                 RequiredNonzeros, context);
              Set<MultipleDeriv> CProd = product(J, K[k], 
                                                 ConstantNonzeros, context);
              SUNDANCE_VERB_EXTREME(tab3 << "CProd = " << CProd);
              SUNDANCE_VERB_EXTREME(tab3 << "RProd = " << RProd);
              Set<MultipleDeriv> X = CProd.setUnion(RM.setDifference(RProd));
              if (!started) 
                {
                  rtn = X;
                  started = true;
                }
              else 
                {
                  rtn = rtn.intersection(X);
                }
            }
        }

      const Set<MultiSet<int> >& QV = findQ_V(i, context);
      SUNDANCE_VERB_EXTREME(tab1 << "finding R\\RProd over QV");      
      SUNDANCE_VERB_EXTREME(tab1 << "QV = " << QV);

      for (Set<MultiSet<int> >::const_iterator j=QV.begin(); j!=QV.end(); j++)
        {
          Tabs tab2;
          Array<int> J = j->elements();
          const Array<Array<int> >& K = comp[J.size()-1];
          SUNDANCE_VERB_EXTREME(tab2 << "J=" << J);
          SUNDANCE_VERB_EXTREME(tab2 << "K=" << K);

          for (unsigned int k=0; k<K.size(); k++)
            {
              Tabs tab3;
              Set<MultipleDeriv> RProd = product(J, K[k], 
                                                 RequiredNonzeros, context);
              Set<MultipleDeriv> X = RM.setDifference(RProd);
              if (!started) 
                {
                  rtn = X;
                  started = true;
                }
              else 
                {
                  rtn = rtn.intersection(X);
                }
              SUNDANCE_VERB_EXTREME(tab3 << "RProd = " << RProd);
            }
        }
    }

  SUNDANCE_VERB_HIGH(tab0 << "C[" << order << "]=" << rtn);
  SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindC for "
                     << toString());
  return rtn;
}



Set<MultipleDeriv> 
ExprWithChildren::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  Set<MultipleDeriv> rtn;
  SUNDANCE_VERB_HIGH(tab0 << "EWC::internalFindW() for " << toString());  
  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
    {
      /* If I am an arbitrary nonlinear expression, I cannot be known to
       * be zero regardless of the state of my arguments. Return the
       * zeroth-order derivative. */
      if (!(isLinear() || isProduct())) 
        {
          rtn.put(MultipleDeriv());
          SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
          SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindW for "
                             << toString());
          return rtn;
        }

      /* At this point, I've dealt with arbitrary nonlinear exprs so 
       * I know I'm either a product or a linear expr */

      const Set<MultiSet<int> >& Q = findQ_W(0, context);

      /* If there are no nonzero terms, a linear combination or a product
       * will be zero. Return the empty set. */
      if (Q.size()==0)
        {
          SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
          SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindW for "
                             << toString());
          return rtn;
        }
      
      /* if I'm a linear combination and any term is nonzero, I am nonzero */
      if (isLinear())
        {
          rtn.put(MultipleDeriv());      
          SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
          SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindW for "
                             << toString());    
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
      SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
      SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindW for "
                         << toString());
      return rtn;
    }


  /* now do arbitrary order derivatives with the multivariable chain rule*/
  Array<Array<Array<int> > > comp = compositions(order);
  for (int i=1; i<=order; i++) 
    {


      const Set<MultiSet<int> >& QW = findQ_W(i, context);
      
      for (Set<MultiSet<int> >::const_iterator j=QW.begin(); j!=QW.end(); j++)
        {
          Array<int> J = j->elements();
          const Array<Array<int> >& K = comp[J.size()-1];

          for (unsigned int k=0; k<K.size(); k++)
            {
              Set<MultipleDeriv> WProd = product(J, K[k], AllNonzeros, context);
              rtn.merge(WProd);
            }
        }
    }

  SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
  SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalFindW for "
                     << toString());
  return rtn;
}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_W(int order, 
                          const EvalContext& context) const
{
  if (!contextToQWMap_[order].containsKey(context))
    {
      contextToQWMap_[order].put(context, internalFindQ_W(order, context));
    }
  return contextToQWMap_[order].get(context);
}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_V(int order, 
                          const EvalContext& context) const
{
  if (!contextToQVMap_[order].containsKey(context))
    {
      contextToQVMap_[order].put(context, internalFindQ_V(order, context));
    }
  return contextToQVMap_[order].get(context);
}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_C(int order, 
                          const EvalContext& context) const
{
  if (!contextToQCMap_[order].containsKey(context))
    {
      contextToQCMap_[order].put(context, internalFindQ_C(order, context));
    }
  return contextToQCMap_[order].get(context);
}


const Set<MultiSet<int> >& ExprWithChildren::getI_N() const
{
  if (!cachedI_N().containsKey(numChildren()))
    {
      Set<MultiSet<int> > x;
      for (int i=0; i<numChildren(); i++)
        {
          x.put(makeMultiSet<int>(i));
        }
      cachedI_N().put(numChildren(), x);
    }
  return cachedI_N().get(numChildren());
}

Set<MultiSet<int> > ExprWithChildren::indexSetProduct(const Set<MultiSet<int> >& a,
                                                      const Set<MultiSet<int> >& b) const
{
  Set<MultiSet<int> > rtn;
  for (Set<MultiSet<int> >::const_iterator i=a.begin(); i!=a.end(); i++)
    {
      for (Set<MultiSet<int> >::const_iterator j=b.begin(); j!=b.end(); j++)
        {
          MultiSet<int> ab = (*i).merge(*j);
          rtn.put(ab);
        }
    }
  return rtn;
}


Set<MultiSet<int> > ExprWithChildren::internalFindQ_V(int order, 
                                                      const EvalContext& context) const
{
  Set<MultiSet<int> > rtn;

  if (!isLinear())
    {
      bool isVar = false;
      for (int i=0; i<numChildren(); i++)
        {
          if (childIsRequired(i,order,context) && evaluatableChild(i)->findV(0, context).size() != 0) 
            {
              isVar=true;
              break;
            }
        }
      if (isVar) rtn = findQ_V(order, context); 
    }
  return rtn;
}

Set<MultiSet<int> > ExprWithChildren::internalFindQ_C(int order, 
                                                      const EvalContext& context) const
{
  if (isLinear()) return findQ_W(order,context);

  return findQ_W(order,context).setDifference(findQ_V(order, context));
}


Set<MultiSet<int> > ExprWithChildren
::internalFindQ_W(int order, 
                  const EvalContext& context) const
{
  Set<MultiSet<int> > rtn;
  const Set<MultiSet<int> >& I_N = getI_N();

  if (isLinear())
    {
      /* first derivatives of the sum wrt the arguments are 
       * always nonzero */
      if (order==1) return I_N;
      /* zeroth derivatives are nonzero if terms are nonzero */
      if (order==0)
        {
          for (int i=0; i<numChildren(); i++)
            {
              const Set<MultipleDeriv>& W_i = evaluatableChild(i)->findW(0,context);
              if (W_i.size() > 0) rtn.put(makeMultiSet(i));
            }
        }
    }
  else
    {
      rtn = I_N;
      
      for (int i=1; i<order; i++)
        {
          rtn = indexSetProduct(rtn, I_N);
        }
    }
  return rtn;
}

bool ExprWithChildren::childIsRequired(int index, int diffOrder,
                                       const EvalContext& context) const
{
  const Set<MultiSet<int> >& Q = findQ_W(diffOrder, context);
  for (Set<MultiSet<int> >::const_iterator it=Q.begin(); it != Q.end(); it++)
    {
      if (it->contains(index)) return true;
    }
  return true;
}

RefCountPtr<Array<Set<MultipleDeriv> > > ExprWithChildren
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab0;
  RefCountPtr<Array<Set<MultipleDeriv> > > rtn 
    = rcp(new Array<Set<MultipleDeriv> >(RInput.size()));

  SUNDANCE_VERB_HIGH(tab0 << "in internalDetermineR() for " << toString());
  SUNDANCE_VERB_HIGH(tab0 << "RInput = " << RInput);



  for (unsigned int i=0; i<RInput.size(); i++)
    {
      const Set<MultipleDeriv>& Wi = findW(i, context);
      (*rtn)[i] = RInput[i].intersection(Wi);
    }

  int maxOrder = rtn->size()-1;

  const Set<MultiSet<int> >& Q1 = findQ_W(1, context);
  const Set<MultiSet<int> >& Q2 = findQ_W(2, context);
  const Set<MultiSet<int> >& Q3 = findQ_W(3, context);

  SUNDANCE_VERB_EXTREME(tab0 << "Q1 = " << Q1);
  SUNDANCE_VERB_EXTREME(tab0 << "Q2 = " << Q2);
  SUNDANCE_VERB_EXTREME(tab0 << "Q3 = " << Q3);

  for (int i=0; i<numChildren(); i++)
    {
      Tabs tab1;
      MultiSet<int> mi = makeMultiSet(i);
      SUNDANCE_VERB_EXTREME(tab1 << "Q1_i = " << mi );
      TEST_FOR_EXCEPTION(mi.size() != 1, InternalError, "unexpected multiset size");
      int i = *(mi.begin());
      Set<MultipleDeriv> R11;
      Set<MultipleDeriv> R12;
      Set<MultipleDeriv> R13;
      Set<MultipleDeriv> R22;
      Set<MultipleDeriv> R23;
      Set<MultipleDeriv> R33;
      if (maxOrder >= 1) 
        {
          R11 = (*rtn)[1];
          if (maxOrder >=2) R22 = (*rtn)[2];
          if (maxOrder >=3) R33 = (*rtn)[3];
        }
      if (maxOrder >= 2)
        {
          Tabs tab2;
          Set<MultiSet<int> > jSet = setDivision(Q2, i);
          SUNDANCE_VERB_EXTREME(tab2 << "Q2/i = " << jSet);
          for (Set<MultiSet<int> >::const_iterator 
                 j=jSet.begin(); j!=jSet.end(); j++)
            {
              Tabs tab3;
              TEST_FOR_EXCEPTION(j->size()!=1, InternalError, 
                                 "unexpected set size");
              int jIndex = *(j->begin());
              SUNDANCE_VERB_EXTREME( tab3 << "j=" << jIndex );
              const Set<MultipleDeriv>& W1j = evaluatableChild(jIndex)->findW(1, context);
              Set<MultipleDeriv> ROverW = setDivision((*rtn)[2], W1j);
              R12.merge(ROverW);
              SUNDANCE_VERB_EXTREME( tab3 << "R2=" << (*rtn)[2] );
              SUNDANCE_VERB_EXTREME( tab3 << "W1(j)=" << W1j );
              SUNDANCE_VERB_EXTREME( tab3 << "R2/W1(j)=" << ROverW );

              if (maxOrder>=3)
                {
                  const Set<MultipleDeriv>& W2j 
                    = evaluatableChild(jIndex)->findW(2, context);             
                  R13.merge(setDivision((*rtn)[3], W2j));
                  R23.merge(setDivision((*rtn)[3], W1j));
                }
            }
        }
      if (maxOrder >= 3)
        {
          Set<MultiSet<int> > jkSet = setDivision(Q3, i);
          for (Set<MultiSet<int> >::const_iterator 
                 jk=jkSet.begin(); jk!=jkSet.end(); jk++)
            {
              TEST_FOR_EXCEPTION(jk->size()!=2, InternalError, 
                                 "unexpected set size");
              Array<int> jka = jk->elements();
              int j = jka[0];
              int k = jka[1];
              const Set<MultipleDeriv>& W1j = evaluatableChild(j)->findW(1, context);
              const Set<MultipleDeriv>& W1k = evaluatableChild(k)->findW(1, context);
              R13.merge(setDivision((*rtn)[3], setProduct(W1j, W1k)));
            }
        }
      SUNDANCE_VERB_EXTREME( tab1 << "R11 = " << R11 );
      SUNDANCE_VERB_EXTREME( tab1 << "R12 = " << R12 );
      SUNDANCE_VERB_EXTREME( tab1 << "R13 = " << R13 );
      SUNDANCE_VERB_EXTREME( tab1 << "R22 = " << R22 );
      SUNDANCE_VERB_EXTREME( tab1 << "R23 = " << R23 );
      SUNDANCE_VERB_EXTREME( tab1 << "R33 = " << R33 );
      Set<MultipleDeriv> R1 = R11;
      R1.merge(R12);
      R1.merge(R13);
      Set<MultipleDeriv> R2 = R22;
      R2.merge(R23);
      Set<MultipleDeriv> R3 = R33;

      Set<MultipleDeriv> R0;
      bool childIsNeeded = (*rtn)[0].size() > 0;
      if (!childIsNeeded && R1.size() > 0) childIsNeeded = childIsRequired(i, 2, context);
      if (!childIsNeeded && R2.size() > 0) childIsNeeded = childIsRequired(i, 3, context);
      if (!childIsNeeded && R3.size() > 0) childIsNeeded = childIsRequired(i, 4, context);
      if (childIsNeeded) R0.put(MultipleDeriv());
      
      Array<Set<MultipleDeriv> > RChild;
      
      RChild.append(R0);
      if (maxOrder >= 1) RChild.append(R1);
      if (maxOrder >= 2) RChild.append(R2);
      if (maxOrder >= 3) RChild.append(R3);
      SUNDANCE_VERB_EXTREME( tab1 << "RChild = " << RChild );
      evaluatableChild(i)->determineR(context, RChild);
    }

  SUNDANCE_VERB_HIGH( tab0 << "R = " << (*rtn) );
  SUNDANCE_VERB_HIGH(tab0 << "done with EWC::internalDetermineR for "
                     << toString());
  
  return rtn;
}




void ExprWithChildren::displayNonzeros(ostream& os, const EvalContext& context) const 
{
  Tabs tabs0;
  os << tabs0 << "Nonzeros of " << toString() << endl;
  os << tabs0 << "Diving into children " << endl;

  for (int i=0; i<numChildren(); i++)
    {
      Tabs tab1;
      os << tab1 << "Child " << i << endl;
      evaluatableChild(i)->displayNonzeros(os, context);
    }

  os << tabs0 << "Printing nonzeros for parent " << toString() << endl;
  const Set<MultipleDeriv>& W = findW(context);
  const Set<MultipleDeriv>& R = findR(context);
  const Set<MultipleDeriv>& C = findC(context);

  
  for (Set<MultipleDeriv>::const_iterator i=W.begin(); i != W.end(); i++)
    {
      Tabs tab1;
      string state = "Variable";
      if (C.contains(*i)) state = "Constant";
      if (!R.contains(*i)) state = "Not Required";
      os << tab1 << std::setw(25) << std::left << i->toString() << ": " << state << endl;
    }
}

