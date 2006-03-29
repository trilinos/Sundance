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

#include "SundanceCoordExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CoordExpr::CoordExpr(int dir, const string& name)
  : LeafExpr(), FuncElementBase(coordName(dir, name), ""), 
    dir_(dir)
{
  setOrderOfDependency(dir, 1);
}

XMLObject CoordExpr::toXML() const 
{
  XMLObject rtn("CoordExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CoordExpr::coordName(int dir, const string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "CoordExpr::coordName direction out of range [0,2]");
      return "error";
    }
}


Set<MultipleDeriv> 
CoordExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  if (order==1) 
    {
      Deriv x = new CoordDeriv(dir_);
      MultipleDeriv md;
      md.put(x);
      rtn.put(md);
    }
  return rtn;
}


void CoordExpr::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                             const Set<MultiSet<int> >& inputActiveFuncIDs,
                             bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for coord func " << toString()
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

  MultiIndex myDirection;
  myDirection[dir_] = 1;
  if (multiIndices.contains(myDirection)) 
    {
      subset->addDeriv(new CoordDeriv(dir_), ConstantDeriv);
     //  if (activeFuncIDs.contains(MultiSet<int>()))
//         {
//           subset->addDeriv(new CoordDeriv(dir_), ConstantDeriv);
//         }
    }
  MultiIndex empty;
  if (multiIndices.contains(empty))
    {
      if (activeFuncIDs.contains(MultiSet<int>()))
        {
          subset->addDeriv(MultipleDeriv(), VectorDeriv);
        }
    }

  

  SUNDANCE_VERB_HIGH(tabs << "coord expr: " + toString() 
                     << ": my sparsity subset is " 
                       << endl << *subset);

  TEST_FOR_EXCEPTION(sparsitySuperset(context).get()==0, InternalError,
                     "null sparsity superset detected in CoordExpr::findNonzeros()");

  SUNDANCE_VERB_HIGH(tabs << "coord expr:  " 
                     + toString() << ": my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}




