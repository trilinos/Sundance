/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

SumOfIntegrals::SumOfIntegrals(const RefCountPtr<CellFilterBase>& domain,
                               const Expr& expr,
                               const RefCountPtr<QuadratureFamilyBase>& quad)
  : ScalarExpr(), domains_(),
    quad_(),
    expr_(),
    cellSetToIndexMap_(),
    quadToIndexMap_()
{
  addTerm(domain, expr, quad, 1);
}


void SumOfIntegrals::addTerm(const RefCountPtr<CellFilterBase>& domainPtr,
                             const Expr& expr,
                             const RefCountPtr<QuadratureFamilyBase>& quadPtr, 
                             int sign)
{
  int d = domains_.length();

  OrderedHandle<CellFilterBase> domain(domainPtr);
  OrderedHandle<QuadratureFamilyBase> quad(quadPtr);

  if (cellSetToIndexMap_.containsKey(domain))
    {
      d = cellSetToIndexMap_.get(domain);
    }
  else
    {
      domains_.append(domain);
      quadToIndexMap_.resize(d+1);
      quad_.resize(d+1);
      expr_.resize(d+1);
      cellSetToIndexMap_.put(domain, d);
    }
  
  int q = quad_[d].length();
  if (quadToIndexMap_[d].containsKey(quad))
    {
      q = quadToIndexMap_[d].get(quad);
      if (sign > 0)
        {
          expr_[d][q] = expr_[d][q] + expr;
        }
      else
        {
          expr_[d][q] = expr_[d][q] - expr;
        }
    }
  else
    {
      quad_[d].append(quad);
      if (sign > 0)
        {
          expr_[d].append(expr);
        }
      else
        {
          expr_[d].append(-expr);
        }
      quadToIndexMap_[d].put(quad, q);
    }
}

void SumOfIntegrals::merge(const SumOfIntegrals* other, int sign) 
{
  cerr << "before merge" << endl;

  {
    Tabs tabs;
    
    cerr << tabs << "me: " << toXML() << endl;
    cerr << tabs << "other: " << other->toXML() << endl;
  }

  for (int d=0; d<other->domains_.size(); d++)
    {
      for (int q=0; q<other->numTerms(d); q++)
        {
          addTerm(other->domains_[d].ptr(), other->expr_[d][q],
                  other->quad_[d][q].ptr(), sign);
        }
    }
  
  cerr << "after merge " << endl;
  {
    Tabs tabs;
    
    cerr << tabs << "me: " << toXML() << endl;
  }
}

void SumOfIntegrals::multiplyByConstant(const SpatiallyConstantExpr* expr) 
{
  double a = expr->value();

  for (int d=0; d<domains_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = a*expr_[d][q];
        }
    }
}

void SumOfIntegrals::changeSign()
{
  for (int d=0; d<domains_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = -expr_[d][q];
        }
    }
}

Set<int> SumOfIntegrals::unksOnDomain(int d) const 
{
  Set<int> rtn;
  for (int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateUnkSet(rtn);
    }
  return rtn;
}

Set<int> SumOfIntegrals::testsOnDomain(int d) const 
{
  Set<int> rtn;
  for (int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateTestSet(rtn);
    }
  return rtn;
}



RefCountPtr<CellFilterBase> SumOfIntegrals::nullDomain() const
{
  for (int d=0; d<domains_.size(); d++)
    {
      if (!domains_[d].ptr()->isNullDomain())
        {
          return domains_[d].ptr()->makeNullDomain();
        }
    }
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "SumOfIntegrals::nullDomain() called on a sum "
                     "of integrals with no non-null domains");

  return RefCountPtr<CellFilterBase>();
}



ostream& SumOfIntegrals::toText(ostream& os, bool paren) const
{
  os << "Sum of Integrals[" << endl;
  for (int d=0; d<domains_.size(); d++)
    {
      for (int t=0; t<quad_[d].size(); t++)
        { 
          os << "Integral[" << endl;
          os << domains_[d].ptr()->toXML() << endl;
          os << "quad rule: " << quad_[d][t].ptr()->toXML() << endl;
          os << "expr: " << expr_[d][t].toString() << endl;
          os << "]" << endl;
        }

    }
  os << "]" << endl;

  return os;
}

ostream& SumOfIntegrals::toLatex(ostream& os, bool paren) const
{
  TEST_FOR_EXCEPTION(true, InternalError, 
                     "SumOfIntegrals::toLatex is undefined");
  return os;
}

XMLObject SumOfIntegrals::toXML() const 
{
  XMLObject rtn("SumOfIntegrals");
  for (int d=0; d<domains_.size(); d++)
    {
      XMLObject child("Integral");
      rtn.addChild(child);
      for (int t=0; t<quad_[d].size(); t++)
        {
          child.addChild(quad_[d][t].ptr()->toXML());
          child.addChild(expr_[d][t].toXML());
        }
    }

  return rtn;
}
