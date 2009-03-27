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
#include "SundanceSpectralExpr.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SumOfIntegrals::SumOfIntegrals(const RefCountPtr<CellFilterStub>& region,
  const Expr& expr,
  const RefCountPtr<QuadratureFamilyStub>& quad,
  const WatchFlag& watch)
  : ScalarExpr(), rqcToExprMap_()
{
  addTerm(region, expr, quad, watch, 1);
}


Expr SumOfIntegrals::filterSpectral(const Expr& expr) const 
{
  const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(expr.ptr().get());
  if (se != 0) return se->getCoeff(0);
  return expr;
}



void SumOfIntegrals::addTerm(const RefCountPtr<CellFilterStub>& regionPtr,
  const Expr& expr,
  const RefCountPtr<QuadratureFamilyStub>& quadPtr, 
  const WatchFlag& watch, int sign)
{
  Expr ex = filterSpectral(expr);

  RegionQuadCombo rqc(regionPtr, quadPtr, watch);

  if (rqcToExprMap_.containsKey(rqc))
  {
    Expr e = rqcToExprMap_.get(rqc);
    rqcToExprMap_.put(rqc, e + sign*expr);
  }
  else
  {
    rqcToExprMap_.put(rqc, sign*expr);
  }
}

void SumOfIntegrals::merge(const SumOfIntegrals* other, int sign) 
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=other->rqcToExprMap_.begin(); i!=other->rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    const Expr& e = i->second;
    addTerm(rqc.domain(), e, rqc.quad(), rqc.watch(), sign);
  }
}

void SumOfIntegrals::multiplyByConstant(const SpatiallyConstantExpr* expr) 
{
  double a = expr->value();
  SundanceUtils::Map<RegionQuadCombo, Expr> newMap;
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    newMap.put(i->first, a*e);
  }
  rqcToExprMap_ = newMap;
}

void SumOfIntegrals::changeSign()
{
  SundanceUtils::Map<RegionQuadCombo, Expr> newMap;
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    newMap.put(i->first, -e);
  }
  rqcToExprMap_ = newMap;
}

Set<int> SumOfIntegrals::funcsOnRegion(const OrderedHandle<CellFilterStub>& d, const Set<int>& funcSet) const 
{
  Set<int> rtn;
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (OrderedHandle<CellFilterStub>(rqc.domain()) != d) continue;
    Expr e = i->second;
    e.ptr()->accumulateFuncSet(rtn, funcSet);
  }
  return rtn;
}


bool SumOfIntegrals::integralHasTestFunctions(const OrderedHandle<CellFilterStub>& d) const 
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (OrderedHandle<CellFilterStub>(rqc.domain()) != d) continue;
    Expr e = i->second;
    if (e.hasTestFunctions()) return true;
  }
  return false;
}



RefCountPtr<CellFilterStub> SumOfIntegrals::nullRegion() const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (!rqc.domain()->isNullRegion())
    {
      return rqc.domain()->makeNullRegion();
    }
  }
  
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "SumOfIntegrals::nullRegion() called on a sum "
                     "of integrals with no non-null regions");

  return RefCountPtr<CellFilterStub>();
}

bool SumOfIntegrals::isIndependentOf(const Expr& u) const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isIndependentOf(u)) return false;
  }
  return true;
}

bool SumOfIntegrals::isLinearForm(const Expr& u) const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isLinearForm(u)) return false;
  }
  return true;
}

bool SumOfIntegrals::isQuadraticForm(const Expr& u) const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isQuadraticForm(u)) return false;
  }
  return true;
}


bool SumOfIntegrals::everyTermHasTestFunctions() const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.everyTermHasTestFunctions()) return false;
  }
  return true;
}


bool SumOfIntegrals::isLinearInTests() const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (!e.isLinearInTests()) return false;
  }
  return true;
}

bool SumOfIntegrals::hasTestFunctions() const
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    Expr e = i->second;
    if (e.hasTestFunctions()) return true;
  }
  return false;
}



ostream& SumOfIntegrals::toText(ostream& os, bool paren) const
{
  os << "Sum of Integrals[" << endl;
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    Expr e = i->second;
    os << "Integral[" << endl;
    os << "rqc=" << rqc.toString() << endl;
    os << "expr=" << e.toString() << endl;
    os << "]" << endl;
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
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    Expr e = i->second;
    XMLObject child("Integral");
    rtn.addChild(child);
    child.addChild(rqc.quad()->toXML());
    child.addChild(rqc.domain()->toXML());
    child.addChild(rqc.watch().toXML());
    child.addChild(e.toXML());
  }

  return rtn;
}


bool SumOfIntegrals::lessThan(const ScalarExpr* other) const
{
  const SumOfIntegrals* f = dynamic_cast<const SumOfIntegrals*>(other);
  TEST_FOR_EXCEPTION(f==0, InternalError, "cast should never fail at this point");
  
  return rqcToExprMap_ < f->rqcToExprMap_;
}

bool SumOfIntegrals::hasWatchedTerm() const 
{
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         i=rqcToExprMap_.begin(); i!=rqcToExprMap_.end(); i++)
  {
    const RegionQuadCombo& rqc = i->first;
    if (rqc.watch().isActive()) return true;
  }
  return false;
}


