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

#include "SundanceLabelCellPredicate.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool LabelCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(dynamic_cast<const LabelCellPredicate*>(other) == 0,
                     std::logic_error,
                     "argument " << other->toXML() 
                     << " to LabelCellPredicate::lessThan() should be "
                     "a LabelCellPredicate pointer.");

  const LabelCellPredicate* op 
    = dynamic_cast<const LabelCellPredicate*>(other);
  
  return (labelIndices_ < op->labelIndices_);
}

void LabelCellPredicate::testBatch(const Array<int>& cellLID,
                                   Array<int>& results) const
{
  mesh().getLabels(cellDim(), cellLID, results);
  for (int i=0; i<cellLID.size(); i++)
    {
      int cellLabel = results[i];
      results[i] = labelIndices_.contains(cellLabel);
    }
}

XMLObject LabelCellPredicate::toXML() const 
{
  XMLObject rtn("LabelCellPredicate");
  rtn.addAttribute("label", labelIndices_.toString());
  return rtn;
}

