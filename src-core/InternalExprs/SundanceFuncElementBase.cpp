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

#include "SundanceFuncElementBase.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

FuncElementBase::FuncElementBase(const string& rootName,
                                 const string& suffix)
	: ScalarExpr(), name_(rootName + suffix), rootName_(rootName),
    suffix_(suffix), id_(nextID()++)
{}

FuncElementBase::FuncElementBase(const string& rootName)
	: ScalarExpr(), name_(rootName), rootName_(rootName),
    suffix_(), id_(nextID()++)
{}

bool FuncElementBase::isInActiveSet(const Set<MultiSet<int> >& activeFuncIDs) const
{
  for (Set<MultiSet<int> >::const_iterator 
         i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      if (i->contains(funcID())) return true;
    }
  return false;
}

void FuncElementBase::accumulateFuncSet(Set<int>& funcIDs,
                                        const Set<int>& activeFuncs) const
{
  if (activeFuncs.contains(funcID())) funcIDs.put(funcID());
}

ostream& FuncElementBase::toText(ostream& os, bool /* paren */) const 
{
	os << name_;
	return os;
}

ostream& FuncElementBase::toLatex(ostream& os, bool /* paren */) const 
{
	os << name_;
	return os;
}



  
