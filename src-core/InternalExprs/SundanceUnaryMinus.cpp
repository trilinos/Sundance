/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryMinus.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

UnaryMinus::UnaryMinus(const RefCountPtr<ScalarExpr>& arg)
  : UnaryExpr(arg)
{;}

ostream& UnaryMinus::toText(ostream& os, bool paren) const 
{
  if (paren) os << "(";
  os << "-" << arg().toString();
  if (paren) os << ")";
  return os;
}

ostream& UnaryMinus::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject UnaryMinus::toXML() const
{
  XMLObject rtn("UnaryMinus");
  rtn.addChild(arg().toXML());
  return rtn;
}