/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNonlinearUnaryOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

NonlinearUnaryOp::NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg,
                                   const RefCountPtr<UnaryFunctor>& op)
  : UnaryExpr(arg), op_(op)
{
  if (isEvaluatable(arg.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          if (evaluatableArg()->orderOfSpatialDependency(d) != 0) 
            {
              setOrderOfDependency(d, -1);
            }
          else
            {
              setOrderOfDependency(d, 0);
            }
        }

      const Set<MultiSet<int> >& argFuncs = evaluatableArg()->funcIDSet();
      typedef Set<MultiSet<int> >::const_iterator iter;

      for (iter i=argFuncs.begin(); i != argFuncs.end(); i++)
        {
          const MultiSet<int>& f1 = *i;
          for (iter j=argFuncs.begin(); j != argFuncs.end(); j++)
            {
              const MultiSet<int>& f2 = *j;
              MultiSet<int> f12 = f1.merge(f2);
              for (iter k=argFuncs.begin(); k != argFuncs.end(); k++)
                {
                  const MultiSet<int>& f3 = *k;
                  
                  if (f1.size()+f2.size()+f3.size() > maxFuncDiffOrder()) 
                    continue;
                  addFuncIDCombo(f12.merge(f3));
                }
            }
        }
    }
}

Set<MultiSet<int> > 
NonlinearUnaryOp::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs) const 
{
  typedef Set<MultiSet<int> >::const_iterator iter;

  Set<MultiSet<int> > rtn;
  for (iter i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& f1 = *i;
      for (iter j=activeFuncIDs.begin(); j != activeFuncIDs.end(); j++)
        {
          const MultiSet<int>& f2 = *j;
          MultiSet<int> f12 = f1.merge(f2);
          for (iter k=activeFuncIDs.begin(); k != activeFuncIDs.end(); k++)
            {
              const MultiSet<int>& f3 = *k;
              
              if (f1.size()+f2.size()+f3.size() > maxFuncDiffOrder()) 
                continue;
              rtn.put(f12.merge(f3));
            }
        }
    }
  return rtn;
}

void NonlinearUnaryOp::findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    const Set<MultiSet<int> >& activeFuncIDs,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for unary nonlinear op " 
                       << toString() << " subject to multiindices "
                       << multiIndices);
  SUNDANCE_VERB_MEDIUM(tabs << "active funcs are " << activeFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  addActiveFuncs(context, activeFuncIDs);
  RefCountPtr<SparsitySubset> subset 
    = sparsitySubset(context, multiIndices, activeFuncIDs);

  Set<MultiSet<int> > childFuncIDs = argActiveFuncs(activeFuncIDs);


  
  SUNDANCE_VERB_MEDIUM(tabs << "arg funcID set is " << childFuncIDs);

  int maxMiOrder = maxOrder(multiIndices);
  int maxDiffOrder = context.topLevelDiffOrder() + maxMiOrder;

  evaluatableArg()->findNonzeros(context, multiIndices,
                                 childFuncIDs,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsitySubset 
    = evaluatableArg()->sparsitySubset(context, multiIndices, childFuncIDs);

  SUNDANCE_VERB_MEDIUM(tabs << "arg sparsity subset is " 
                       << endl << *argSparsitySubset);

  for (int i=0; i<argSparsitySubset->numDerivs(); i++)
    {
      if (argSparsitySubset->deriv(i).order()==0)
        {
          subset->addDeriv(argSparsitySubset->deriv(i), 
                           argSparsitySubset->state(i));
        }
      else
        {
          subset->addDeriv(argSparsitySubset->deriv(i), 
                           VectorDeriv);
        }
    }

  for (int i=0; i<argSparsitySubset->numDerivs(); i++)
    {
      for (int j=0; j<argSparsitySubset->numDerivs(); j++)
        {
          MultipleDeriv product 
            = argSparsitySubset->deriv(i).product(argSparsitySubset->deriv(j));
          if (product.order() > maxDiffOrder) continue;
          if (product.spatialOrder() > maxMiOrder) continue;
          subset->addDeriv(product, VectorDeriv);
        }
    }



  SUNDANCE_VERB_HIGH(tabs << "nonlinear op  " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "nonlinear op " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}


ostream& NonlinearUnaryOp::toText(ostream& os, bool paren) const 
{
  os << op_->name() << "(" << arg().toString() << ")";
  return os;
}

ostream& NonlinearUnaryOp::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject NonlinearUnaryOp::toXML() const
{
  XMLObject rtn("NonlinearUnaryOp");
  rtn.addChild(arg().toXML());
  rtn.addAttribute("op", op_->name());
  return rtn;
}

