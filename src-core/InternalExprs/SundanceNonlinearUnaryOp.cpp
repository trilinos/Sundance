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
          if (evaluatableArg()->orderOfDependency(d) != 0) 
            {
              setOrderOfDependency(d, -1);
            }
          else
            {
              setOrderOfDependency(d, 0);
            }
        }
    }
}

void NonlinearUnaryOp::findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for unary nonlinear op" 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  if (nonzerosAreKnown(context, multiIndices, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  int maxMiOrder = maxOrder(multiIndices);
  int maxDiffOrder = context.topLevelDiffOrder() + maxMiOrder;

  evaluatableArg()->findNonzeros(context, multiIndices,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsitySubset 
        = evaluatableArg()->sparsitySubset(context, multiIndices);

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

  addKnownNonzero(context, multiIndices, regardFuncsAsConstant);
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

