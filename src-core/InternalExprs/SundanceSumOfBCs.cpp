/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumOfBCs.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SumOfBCs::SumOfBCs(const RefCountPtr<CellFilterStub>& region,
                   const Expr& expr,
                   const RefCountPtr<QuadratureFamilyStub>& quad)
  : SumOfIntegrals(region, expr, quad)
{;}





ostream& SumOfBCs::toText(ostream& os, bool paren) const
{
  os << "Sum of BCs[" << endl;
  for (int d=0; d<numRegions(); d++)
    {
      for (int t=0; t<numTerms(d); t++)
        {
          os << "BC[" << endl;
          os << region(d)->toXML() << endl;
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
  for (int d=0; d<numRegions(); d++)
    {
      rtn.addChild(region(d)->toXML());
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
