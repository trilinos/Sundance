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

#include "SundanceCellSet.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"
#include <algorithm>
#include <iterator>

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellSet CellSet::setUnion(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("union", other);
  
  Set<int>& cells = rtn->cells();

  std::set_union(begin(), end(), other.begin(), other.end(), 
                 std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setIntersection(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("intersection", other);
  
  Set<int>& cells = rtn->cells();

  std::set_intersection(begin(), end(), other.begin(), other.end(), 
                        std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setDifference(const CellSet& other) const
{
  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("difference", other);
  
  Set<int>& cells = rtn->cells();

  std::set_difference(begin(), end(), other.begin(), other.end(), 
                      std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}


void CellSet::checkCompatibility(const string& op, const CellSet& other) const 
{
  TEST_FOR_EXCEPTION(meshID() != other.meshID(), RuntimeError,
                     "CellSet::checkCompatibility(): "
                     "incompatible mesh ID numbers in " << op
                     << ". LHS=" << meshID() << " RHS=" << other.meshID());

  TEST_FOR_EXCEPTION(dimension() != other.dimension(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible dimensions in " << op
                     << "LHS has "
                     "dimension=" << dimension() << " but RHS has dimension="
                     << other.dimension());
  
  TEST_FOR_EXCEPTION(cellType() != other.cellType(), RuntimeError,
                     "CellSet::checkCompatibility() incompatible cell types. "
                     " in " << op << " LHS has "
                     "cellType=" << cellType() << " but RHS has cellType="
                     << other.cellType());

  SUNDANCE_OUT(verbosity() > VerbMedium,
               "Set operation: " << op << endl
               << "LHS cells: " << *this << endl
               << "RHS cells: " << other);
               
}



