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

#include "SundanceDeriv.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
 

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


bool Deriv::operator<(const Deriv& other) const 
{
  return ptr()->lessThan(other);
}

bool Deriv::operator==(const Deriv& other) const
{
  if (*this < other || other < *this) return false;
  return true;
}


string Deriv::toString() const 
{
  return ptr()->toString();
}

const FunctionalDeriv* Deriv::funcDeriv() const 
{
  return dynamic_cast<FunctionalDeriv*>(ptr().get());
}

const CoordDeriv* Deriv::coordDeriv() const 
{
  return dynamic_cast<CoordDeriv*>(ptr().get());
}

bool Deriv::isTestFunction() const 
{
  if (!isFunctionalDeriv()) return false;

  const FunctionalDeriv* f = funcDeriv();

  TEST_FOR_EXCEPTION(f==0, RuntimeError,
                     "contents not a functional derivative in Deriv::isTestFunction()");

  const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(f->func());
  return t != 0;
}

bool Deriv::isUnknownFunction() const 
{
  if (!isFunctionalDeriv()) return false;

  const FunctionalDeriv* f = funcDeriv();

  TEST_FOR_EXCEPTION(f==0, RuntimeError,
                     "contents not a functional derivative in Deriv::isTestFunction()");

  const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(f->func());
  return u != 0;
}
