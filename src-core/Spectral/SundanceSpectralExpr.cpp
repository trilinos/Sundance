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

#include "SundanceSpectralExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
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

SpectralBasis SpectralExpr::getSpectralBasis() const
{ 
  return *sbasis_;
}


Expr SpectralExpr::getCoeff(int i) const
{
  return coeffs_[i];
}


ostream& SpectralExpr::toText(ostream& os, bool paren) const
{
  os << "{";
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
