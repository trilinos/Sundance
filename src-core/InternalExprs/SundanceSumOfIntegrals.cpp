/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SumOfIntegrals::SumOfIntegrals(const RefCountPtr<CellFilterStub>& region,
                               const Expr& expr,
                               const RefCountPtr<QuadratureFamilyStub>& quad)
  : ScalarExpr(), regions_(),
    quad_(),
    expr_(),
    cellSetToIndexMap_(),
    quadToIndexMap_()
{
  addTerm(region, expr, quad, 1);
}


void SumOfIntegrals::addTerm(const RefCountPtr<CellFilterStub>& regionPtr,
                             const Expr& expr,
                             const RefCountPtr<QuadratureFamilyStub>& quadPtr, 
                             int sign)
{
  int d = regions_.length();

  OrderedHandle<CellFilterStub> region(regionPtr);
  OrderedHandle<QuadratureFamilyStub> quad(quadPtr);

  if (cellSetToIndexMap_.containsKey(region))
    {
      d = cellSetToIndexMap_.get(region);
    }
  else
    {
      regions_.append(region);
      quadToIndexMap_.resize(d+1);
      quad_.resize(d+1);
      expr_.resize(d+1);
      cellSetToIndexMap_.put(region, d);
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
  for (int d=0; d<other->regions_.size(); d++)
    {
      for (int q=0; q<other->numTerms(d); q++)
        {
          addTerm(other->regions_[d].ptr(), other->expr_[d][q],
                  other->quad_[d][q].ptr(), sign);
        }
    }
}

void SumOfIntegrals::multiplyByConstant(const SpatiallyConstantExpr* expr) 
{
  double a = expr->value();

  for (int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = a*expr_[d][q];
        }
    }
}

void SumOfIntegrals::changeSign()
{
  for (int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = -expr_[d][q];
        }
    }
}

Set<int> SumOfIntegrals::unksOnRegion(int d) const 
{
  Set<int> rtn;
  for (int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateUnkSet(rtn);
    }
  return rtn;
}

Set<int> SumOfIntegrals::testsOnRegion(int d) const 
{
  Set<int> rtn;
  for (int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateTestSet(rtn);
    }
  return rtn;
}



RefCountPtr<CellFilterStub> SumOfIntegrals::nullRegion() const
{
  for (int d=0; d<regions_.size(); d++)
    {
      if (!regions_[d].ptr()->isNullRegion())
        {
          return regions_[d].ptr()->makeNullRegion();
        }
    }
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "SumOfIntegrals::nullRegion() called on a sum "
                     "of integrals with no non-null regions");

  return RefCountPtr<CellFilterStub>();
}



ostream& SumOfIntegrals::toText(ostream& os, bool paren) const
{
  os << "Sum of Integrals[" << endl;
  for (int d=0; d<regions_.size(); d++)
    {
      for (int t=0; t<quad_[d].size(); t++)
        { 
          os << "Integral[" << endl;
          os << regions_[d].ptr()->toXML() << endl;
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
  for (int d=0; d<regions_.size(); d++)
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
