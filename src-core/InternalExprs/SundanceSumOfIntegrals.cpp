/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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
  for (unsigned int d=0; d<other->regions_.size(); d++)
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

  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = a*expr_[d][q];
        }
    }
}

void SumOfIntegrals::changeSign()
{
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (int q=0; q<numTerms(d); q++)
        {
          expr_[d][q] = -expr_[d][q];
        }
    }
}

Set<int> SumOfIntegrals::funcsOnRegion(int d, const Set<int>& funcSet) const 
{
  Set<int> rtn;
  for (unsigned int t=0; t<expr_[d].size(); t++)
    { 
      expr_[d][t].ptr()->accumulateFuncSet(rtn, funcSet);
    }
  return rtn;
}


bool SumOfIntegrals::integralHasTestFunctions(int d) const 
{
  for (unsigned int t=0; t<expr_[d].size(); t++)
    { 
      if (expr_[d][t].ptr()->hasTestFunctions()) return true;
    }
  return false;
}



RefCountPtr<CellFilterStub> SumOfIntegrals::nullRegion() const
{
  for (unsigned int d=0; d<regions_.size(); d++)
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
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      for (unsigned int t=0; t<quad_[d].size(); t++)
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
  for (unsigned int d=0; d<regions_.size(); d++)
    {
      XMLObject child("Integral");
      rtn.addChild(child);
      for (unsigned int t=0; t<quad_[d].size(); t++)
        {
          child.addChild(quad_[d][t].ptr()->toXML());
          child.addChild(expr_[d][t].toXML());
        }
    }

  return rtn;
}
