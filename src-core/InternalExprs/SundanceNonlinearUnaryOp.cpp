/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNonlinearUnaryOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

NonlinearUnaryOp::NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg,
                                     const RefCountPtr<UnaryFunctor>& op)
  : UnaryExpr(arg), op_(op)
{;}

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

bool NonlinearUnaryOp::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  MultipleDeriv::const_iterator iter;
  
  for (iter=d.begin(); iter!=d.end(); iter++)
    {
      MultipleDeriv single;
      single.put(*iter);
      if (!evaluatableArg()->hasNonzeroDeriv(single)) return false;
    }
  return true;
}
