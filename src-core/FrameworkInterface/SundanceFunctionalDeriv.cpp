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

#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceOrderedTuple.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

FunctionalDeriv::FunctionalDeriv(const FuncElementBase* func,
                                 const MultiIndex& mi)
  : DerivBase(), func_(func), mi_(mi)
{;}

bool FunctionalDeriv::lessThan(const Deriv& other) const
{
  /* First compare types: spatial derivs are ranked lower than 
   * functional derivs, so if the other guy is a spatial deriv, I'm higher */
  if (other.isCoordDeriv()) return false;

  /* We can now safely cast to a functional deriv and compare contents */

  const FunctionalDeriv* f = other.funcDeriv();

  if (funcID() < f->funcID()) return true;
  if (funcID() > f->funcID()) return false;
  if (multiIndex() < f->multiIndex()) return true;
  return false;
}

Deriv FunctionalDeriv::derivWrtMultiIndex(const MultiIndex& mi) const
{
  return new FunctionalDeriv(func_, mi_ + mi);
}

string FunctionalDeriv::toString() const 
{
  if (mi_.order()==0)
    {
      return func_->name();
    }
  return "D[" + func_->name() + ", " + mi_.coordForm() + "]";
}


