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

#include "SundanceDefs.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceComplexExpr.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Array.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore;
using namespace Teuchos;
using namespace TSFExtended;


SpectralExpr::SpectralExpr(const SpectralBasis& sbasis, const Array<Expr>& coeffs)
  : ExprBase(), 
    coeffs_(),
    sbasis_()

{
  int nterms = sbasis.nterms();
  coeffs_.resize(nterms);
  for (int i=0; i<nterms; i++)
    coeffs_[i] = coeffs[i];
  sbasis_ = rcp(new SpectralBasis(sbasis));
}


SpectralExpr::SpectralExpr(const SpectralBasis& sbasis, const Expr& coeffs)
  : ExprBase(), 
    coeffs_(),
    sbasis_()

{
  int nterms = sbasis.nterms();
  coeffs_.resize(nterms);
  for (int i=0; i<nterms; i++)
    coeffs_[i] = coeffs[i];
  sbasis_ = rcp(new SpectralBasis(sbasis));
}

void SpectralExpr::accumulateFuncSet(Set<int>& funcDofIDs, 
  const Set<int>& activeSet) const
{
  for (unsigned int i=0; i<coeffs_.size(); i++)
  {
    dynamic_cast<const ScalarExpr*>(coeffs_[i].ptr().get())->accumulateFuncSet(funcDofIDs, activeSet);
  }
}

SpectralBasis SpectralExpr::getSpectralBasis() const
{ 
  return *sbasis_;
}


Expr SpectralExpr::getCoeff(int i) const
{
  return coeffs_[i];
}

Expr SpectralExpr::spectralDotProduct(const SpectralExpr* other) const
{
  Expr rtn = coeffs_[0] * other->coeffs_[0];
  for (unsigned int i=1; i<coeffs_.size(); i++) rtn = rtn + coeffs_[i]*other->coeffs_[i];
  return rtn;
}

bool SpectralExpr::hasTestFunctions() const
{
  bool rtn = coeffs_[0].ptr()->hasTestFunctions();
  for (unsigned int i=1; i<coeffs_.size(); i++)
    {
      TEST_FOR_EXCEPTION(coeffs_[i].ptr()->hasTestFunctions() != rtn,
                         InternalError,
                         "expr " << toString() << " has a mix of test and "
                         "non-test coefficients");
    }
  return rtn;
}

bool SpectralExpr::hasHungryDiffOp() const
{
  for (unsigned int i=0; i<coeffs_.size(); i++)
    {
      Expr re = coeffs_[i].real();
      Expr im = coeffs_[i].imag();
      const ScalarExpr* sr = dynamic_cast<const ScalarExpr*>(re.ptr().get());
      const ScalarExpr* si = dynamic_cast<const ScalarExpr*>(im.ptr().get());
      TEST_FOR_EXCEPTION(sr == 0 || si == 0, InternalError,
                         "spectral expr " << toString() << " contains a "
                         "non-scalar coefficient");
      if (sr->isHungryDiffOp() || si->isHungryDiffOp()) return true;
    }
  return false;
}


ostream& SpectralExpr::toText(ostream& os, bool paren) const
{
  os << "SpectralExpr{";
  for (unsigned int i=0; i<coeffs_.size(); i++)
    {
      coeffs_[i].ptr()->toText(os, paren);
      if (i < coeffs_.size()-1) os << ", ";
    }
  os << "}";
  return os;
}

ostream& SpectralExpr::toLatex(ostream& os, bool paren) const
{
  os << "\\{";
  for (unsigned int i=0; i<coeffs_.size(); i++)
    {
      coeffs_[i].ptr()->toLatex(os, paren);
      if (i < coeffs_.size()-1) os << ", ";
    }
  os << "\\}";
  return os;
}

XMLObject SpectralExpr::toXML() const 
{
  XMLObject rtn("SpectralExpr");
  for (int i=0; i<coeffs_.length(); i++)
    {
      rtn.addChild(coeffs_[i].toXML());
    }
  return rtn;
}


namespace SundanceCore
{
  /** */
  Expr getSpectralCoeff(int i, const Expr& e)
  {
    const SpectralExpr* s 
      = dynamic_cast<const SpectralExpr*>(e[0].ptr().get());
    TEST_FOR_EXCEPTION(s!=0, RuntimeError,
                       "getSpectralCoeff() called on non-spectral expr "
                       << e.toString());
    return s->getCoeff(i);
  }

  /** */
  SpectralBasis getSpectralBasis(const Expr& e) 
  {
    const SpectralExpr* s 
      = dynamic_cast<const SpectralExpr*>(e[0].ptr().get());
    TEST_FOR_EXCEPTION(s!=0, RuntimeError,
                       "getSpectralBasis() called on non-spectral expr "
                       << e.toString());
    return s->getSpectralBasis();
  }
}
