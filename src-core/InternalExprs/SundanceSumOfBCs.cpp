/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumOfBCs.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

SumOfBCs::SumOfBCs(const RefCountPtr<CellFilterBase>& domain,
                   const Expr& expr,
                   const RefCountPtr<QuadratureFamilyBase>& quad)
  : SumOfIntegrals(domain, expr, quad)
{;}





ostream& SumOfBCs::toText(ostream& os, bool paren) const
{
  os << "Sum of BCs[" << endl;
  for (int d=0; d<numDomains(); d++)
    {
      for (int t=0; t<numTerms(d); t++)
        {
          os << "BC[" << endl;
          os << domain(d)->toXML() << endl;
          os << "quad rule: " << quad(d,t)->toXML() << endl;
          os << "expr: " << expr(d,t).toString() << endl;
          os << "]" << endl;
        }
    }
  os << "]" << endl;

  return os;
}

ostream& SumOfBCs::toLatex(ostream& os, bool paren) const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "SumOfIntegrals::toLatex is undefined");
  return os;
}

XMLObject SumOfBCs::toXML() const 
{
  XMLObject rtn("SumOfBCs");
  for (int d=0; d<numDomains(); d++)
    {
      rtn.addChild(domain(d)->toXML());
      XMLObject child("BC");
      rtn.addChild(child);
      for (int t=0; t<numTerms(d); t++)
        {
          child.addChild(quad(d,t)->toXML());
          child.addChild(expr(d,t).toXML());
        }
    }

  return rtn;
}
