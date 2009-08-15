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

#include "SundanceCellFilterBase.hpp"
#include "SundanceOut.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;


CellFilterBase::CellFilterBase()
  : CellFilterStub(), cellSetCache_(), subsets_(), disjoints_(), name_()
{;}

CellSet CellFilterBase::getCells(const Mesh& mesh) const
{
  if (cellSetCache_.ptr().get()==0)
    {
      cellSetCache_ = internalGetCells(mesh);
    }
  return cellSetCache_;
}

void CellFilterBase::registerSubset(const CellFilter& sub) const
{
  subsets_.put(sub);
  for (Set<CellFilter>::const_iterator 
         i=sub.knownSubsets().begin(); i!=sub.knownSubsets().end(); i++)
    {
      subsets_.put(*i);
    }
}

void CellFilterBase::registerLabeledSubset(int label, const CellFilter& sub) const
{
  labeledSubsets_.put(label, sub);
  for (SundanceUtils::Map<int, CellFilter>::const_iterator 
         iter=labeledSubsets_.begin(); iter != labeledSubsets_.end(); iter++)
    {
      if (iter->first == label) continue;
      sub.cfbPtr()->registerDisjoint(iter->second);
      iter->second.cfbPtr()->registerDisjoint(sub);
    }
}

void CellFilterBase::registerDisjoint(const CellFilter& sub) const
{
  disjoints_.put(sub);
  for (Set<CellFilter>::const_iterator 
         i=sub.knownSubsets().begin(); i!=sub.knownSubsets().end(); i++)
    {
      disjoints_.put(*i);
    }
}

