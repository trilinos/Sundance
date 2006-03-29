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

#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellDiameterExpr::CellDiameterExpr(const string& name)
  : LeafExpr(), name_(name)
{}

XMLObject CellDiameterExpr::toXML() const 
{
  XMLObject rtn("CellDiameterExpr");
  rtn.addAttribute("name", name_);
  return rtn;
}



Set<MultipleDeriv> 
CellDiameterExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  
  if (order==0) rtn.put(MultipleDeriv());
  
  return rtn;
}



void CellDiameterExpr::findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    const Set<MultiSet<int> >& inputActiveFuncIDs,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for cell diameter " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());

  Set<MultiSet<int> > activeFuncIDs = filterActiveFuncs(inputActiveFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs, false);

  MultiIndex empty;
  if (multiIndices.contains(empty))
    {
      subset->addDeriv(MultipleDeriv(), VectorDeriv);
    }

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}

ostream& CellDiameterExpr::toText(ostream& os, bool paren) const
{
  os << name();
  return os;
}


ostream& CellDiameterExpr::toLatex(ostream& os, bool paren) const
{
  os << name();
  return os;
}


