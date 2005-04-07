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

#include "SundanceCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool CellFilter::isNullCellFilter() const 
{
  return dynamic_cast<const NullCellFilterStub*>(ptr().get()) != 0;
}

CellSet CellFilter::getCells(const Mesh& mesh) const
{
  if (isNullCellFilter())
    {
      return new ExplicitCellSet(mesh, -1, 
                                 NullCell);
    }
  return cfbPtr()->getCells(mesh);
}



int CellFilter::dimension(const Mesh& mesh) const
{
  if (isNullCellFilter())
    {
      return -1;
    }
  return cfbPtr()->dimension(mesh);
}



CellFilter CellFilter::operator+(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Union);
}



CellFilter CellFilter::operator-(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Difference);
}



CellFilter CellFilter::intersection(const CellFilter& other) const 
{
  return new BinaryCellFilter(*this, other, BinaryCellFilter::Intersection);
}



CellFilter CellFilter::labeledSubset(int label) const
{
  CellPredicate pred = new LabelCellPredicate(label);
  return new SubsetCellFilter(*this, pred);
}


CellFilter CellFilter::subset(const CellPredicate& pred) const
{
  return new SubsetCellFilter(*this, pred);
}


CellFilter CellFilter::subset(const RefCountPtr<CellPredicateFunctorBase>& test) const
{
  return new SubsetCellFilter(*this, CellPredicate(test));
}



XMLObject CellFilter::toXML() const 
{
  return ptr()->toXML();
}

const CellFilterBase* CellFilter::cfbPtr() const
{
  const CellFilterBase* rtn = dynamic_cast<CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::cfbPtr() cast failed");
  return rtn;
}

