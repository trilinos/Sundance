/* @HEADER@ */
/* @HEADER@ */

#include "SundanceListExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


ListExpr::ListExpr()
  : ExprBase(), elements_()
{;}

ListExpr::ListExpr(const Array<Expr>& elements)
  : ExprBase(), elements_(elements)
{;}

void ListExpr::append(const Expr& expr)
{
  elements_.append(expr);
}

Expr ListExpr::flatten() const 
{
  Expr rtn = new ListExpr();

  for (int i=0; i<size(); i++)
    {
      Expr e = element(i).flatten();
      for (int j=0; j<e.size(); j++)
        {
          rtn.append(e[j]);
        }
    }

  return rtn;
}

Expr ListExpr::join(const Expr& other) const 
{
  Expr rtn = new ListExpr(elements_);
  
  for (int i=0; i<other.size(); i++)
    {
      rtn.append(other[i]);
    }

  return rtn;
}

int ListExpr::size() const
{
  return elements_.size();
}

int ListExpr::totalSize() const 
{
  int rtn = 0;

  for (int i=0; i<size(); i++)
    {
      rtn += elements_[i].totalSize();
    }

  return rtn;
}

ostream& ListExpr::toText(ostream& os, bool paren) const
{
  os << "{";
  for (int i=0; i<elements_.size(); i++)
    {
      elements_[i].ptr()->toText(os, paren);
      if (i < elements_.size()-1) os << ", ";
    }
  os << "}";
  return os;
}

ostream& ListExpr::toLatex(ostream& os, bool paren) const
{
  os << "\\{";
  for (int i=0; i<elements_.size(); i++)
    {
      elements_[i].ptr()->toLatex(os, paren);
      if (i < elements_.size()-1) os << ", ";
    }
  os << "\\}";
  return os;
}

XMLObject ListExpr::toXML() const 
{
  XMLObject rtn("ListExpr");
  for (int i=0; i<elements_.length(); i++)
    {
      rtn.addChild(elements_[i].toXML());
    }
  return rtn;
}


