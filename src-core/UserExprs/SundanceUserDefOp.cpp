/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUserDefOp.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

UserDefOp::UserDefOp(const Expr& arg,
                     const RefCountPtr<UserDefFunctor>& op)
  : ExprWithChildren(getScalarArgs(arg)), op_(op)
{
  for (int i=0; i<numChildren(); i++)
    {
      if (isEvaluatable(evaluatableChild(i)))
        {
          for (int d=0; d<MultiIndex::maxDim(); d++) 
            {
              if (evaluatableChild(i)->orderOfDependency(d) != 0) 
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
}

ostream& UserDefOp::toText(ostream& os, bool paren) const 
{
  os << op_->name() << "(";
  for (int i=0; i<numChildren(); i++)
    {
      os << child(i).toString();
      if (i < numChildren()-1) os << ",";
    }
  os << ")";
  return os;
}

ostream& UserDefOp::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject UserDefOp::toXML() const
{
  XMLObject rtn("UserDefOp");
  XMLObject args("Arguments");
  for (int i=0; i<numChildren(); i++)
    {
      args.addChild(child(i).toXML());
    }
  rtn.addChild(args);
  rtn.addAttribute("op", op_->name());
  return rtn;
}

void UserDefOp::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                             bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for user-defined "
                       "nonlinear op" << toString() 
                       << " subject to multiindices " << multiIndices); 



  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  Set<MultiSet<int> > childFuncIDs = findChildFuncIDSet(activeFuncIDs,
                                                        allFuncIDs);

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  for (int i=0; i<numChildren(); i++)
    {
      evaluatableChild(i)->findNonzeros(context, multiIndices,
                                        childFuncIDs,
                                        allFuncIDs,
                                        regardFuncsAsConstant);

      RefCountPtr<SparsitySubset> argSparsitySubset 
        = evaluatableChild(i)->sparsitySubset(context, multiIndices);

      for (int j=0; j<argSparsitySubset->numDerivs(); j++)
        {
          subset->addDeriv(argSparsitySubset->deriv(j), VectorDeriv);
        }
    }

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant);
}


Array<RefCountPtr<ScalarExpr> > UserDefOp::getScalarArgs(const Expr& args)
{
  Expr fargs = args.flatten();
  Array<RefCountPtr<ScalarExpr> > sargs(fargs.size());
  
  for (int i=0; i<fargs.size(); i++)
    {
      sargs[i] = rcp_dynamic_cast<ScalarExpr>(args[i].ptr());
    }
  return sargs;
}
