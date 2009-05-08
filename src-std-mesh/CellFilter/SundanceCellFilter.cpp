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
#include "SundanceCellFilterBase.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceNullCellFilter.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace Teuchos;

bool CellFilter::isNullCellFilter() const 
{
  return dynamic_cast<const NullCellFilterStub*>(ptr().get()) != 0;
}

bool CellFilter::isNull() const
{
  return ptr().get() == 0 || isNullCellFilter();
}

void CellFilter::setName(const string& name)
{
  nonConstCfbPtr()->setName(name);
}

CellSet CellFilter::getCells(const Mesh& mesh) const
{
  if (isNull() || isNullCellFilter())
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
  if (isNull())
    {
      return other;
    }
  else if (other.isNull())
    {
      return *this;
    }
  else
    {
      CellFilter rtn 
        = new BinaryCellFilter(*this, other, BinaryCellFilter::Union);
      rtn.cfbPtr()->registerSubset(*this);
      rtn.cfbPtr()->registerSubset(other);
      return rtn;
    }
}



CellFilter CellFilter::operator-(const CellFilter& other) const 
{
  if (other.isNull())
    {
      return *this;
    }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
    {
      return *this;
    }
  else if (isKnownSubsetOf(other))
    {
      CellFilter rtn;
      return rtn;
    }
  else if (*this == other)
    {
      CellFilter rtn;
      return rtn;
    }
  else
    {
      CellFilter rtn 
        = new BinaryCellFilter(*this, other, BinaryCellFilter::Difference);
      rtn.cfbPtr()->registerDisjoint(other);
      this->cfbPtr()->registerSubset(rtn);
      return rtn;
    }
}



CellFilter CellFilter::intersection(const CellFilter& other) const 
{
  if (isNull() || other.isNull())
    {
      CellFilter rtn;
      return rtn;
    }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
    {
      CellFilter rtn;
      return rtn;
    }
  else if (isKnownSubsetOf(other))
    {
      return *this;
    }
  else if (other.isKnownSubsetOf(*this))
    {
      return other;
    }
  else if (*this==other)
    {
      return *this;
    }
  else
    {
      CellFilter rtn 
        = new BinaryCellFilter(*this, other, BinaryCellFilter::Intersection);
      this->cfbPtr()->registerSubset(rtn);
      other.cfbPtr()->registerSubset(rtn);
      return rtn;
    }
}



CellFilter CellFilter::labeledSubset(int label) const
{
  CellPredicate pred = new LabelCellPredicate(label);
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  cfbPtr()->registerLabeledSubset(label, rtn);
  cfbPtr()->registerSubset(rtn);
  return rtn;
}

CellFilter CellFilter::subset(const CellPredicate& pred) const
{
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  cfbPtr()->registerSubset(rtn);
  return rtn;
}


CellFilter CellFilter::subset(const RefCountPtr<CellPredicateFunctorBase>& test) const
{
  CellFilter rtn = new SubsetCellFilter(*this, CellPredicate(test));
  cfbPtr()->registerSubset(rtn);
  return rtn;
}

const Set<CellFilter>& CellFilter::knownSubsets() const
{
  return cfbPtr()->knownSubsets();
}

const Set<CellFilter>& CellFilter::knownDisjoints() const
{
  return cfbPtr()->knownDisjoints();
}

bool CellFilter::isKnownSubsetOf(const CellFilter& other) const
{
  if (other.knownSubsets().contains(*this)) return true;
  return false;
}

bool CellFilter::isKnownDisjointWith(const CellFilter& other) const
{
  if (other.knownDisjoints().contains(*this)) return true;
  if (this->knownDisjoints().contains(other)) return true;

  return false;
}

bool CellFilter::isSubsetOf(const CellFilter& other,
                            const Mesh& mesh) const
{
  if (isKnownSubsetOf(other)) 
    {
      return true;
    }
  else
    {
      CellSet myCells = getCells(mesh);
      CellSet yourCells = other.getCells(mesh);
      CellSet inter = myCells.setIntersection(yourCells);
      if (inter.begin() == inter.end()) return false;
      CellSet diff = myCells.setDifference(inter);
      return (diff.begin() == diff.end());
    }
}



bool CellFilter::operator==(const CellFilter& other) const
{
  if (*this < other) return false;
  if (other < *this) return false;
  return true;
}

bool CellFilter::operator!=(const CellFilter& other) const
{
  return !( *this == other );
}

XMLObject CellFilter::toXML() const 
{
  return ptr()->toXML();
}

string CellFilter::toString() const 
{
  return cfbPtr()->toString();
}

const CellFilterBase* CellFilter::cfbPtr() const
{
  const CellFilterBase* rtn = dynamic_cast<const CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::cfbPtr() cast failed");
  return rtn;
}

CellFilterBase* CellFilter::nonConstCfbPtr()
{
  CellFilterBase* rtn = dynamic_cast<CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::nonConstCfbPtr() cast failed");
  return rtn;
}
