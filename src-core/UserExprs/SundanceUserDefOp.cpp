/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUserDefOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

UserDefOp::UserDefOp(const Expr& arg,
                     const RefCountPtr<UserDefFunctor>& op)
  : ExprWithChildren(getScalarArgs(arg)), op_(op)
{;}

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

bool UserDefOp::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  TimeMonitor t(nonzeroDerivCheckTimer());
  hasNonzeroDerivCalls()++;

  if (derivHasBeenCached(d))
    {
      nonzeroDerivCacheHits()++;
      return getCachedDerivNonzeroness(d);
    }
  TimeMonitor t2(uncachedNonzeroDerivCheckTimer());

  MultipleDeriv::const_iterator iter;
  
 

  bool gotit = false;
  for (int i=0; i<numChildren(); i++)
    {
      if (evaluatableChild(i)->hasNonzeroDeriv(d))
        {
          addDerivToCache(d, true);
          gotit = true;
        }
    }
  if (gotit) return true;

  addDerivToCache(d, false);
  return false;
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
