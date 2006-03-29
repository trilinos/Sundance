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

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

DiscreteFuncElement
::DiscreteFuncElement(const RefCountPtr<DiscreteFuncDataStub>& data,
                      const string& name,
                      const string& suffix,
                      int myIndex)
	: LeafExpr(), 
    FuncElementBase(name, suffix),
    commonData_(data),
    myIndex_(myIndex)
{}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) 
    {
      rtn.put(MultipleDeriv());
    }
  if (order==1)
    {
      Deriv x = new CoordDeriv(0);
      Deriv y = new CoordDeriv(1);
      Deriv z = new CoordDeriv(2);
      MultipleDeriv mx;
      MultipleDeriv my;
      MultipleDeriv mz;
      mx.put(x);
      my.put(y);
      mz.put(z);
      rtn.put(mx);
      rtn.put(my);
      rtn.put(mz);
    }

  return rtn;
}

void DiscreteFuncElement::findNonzeros(const EvalContext& context,
                                       const Set<MultiIndex>& multiIndices,
                                       const Set<MultiSet<int> >& inputActiveFuncIDs,
                                       bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "DiscreteFunc " << toString() << " finding nonzeros " 
                       " subject to multiindices "
                       << multiIndices);

  Set<MultiSet<int> > activeFuncIDs = filterActiveFuncs(inputActiveFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs, false);
  
  for (Set<MultiIndex>::const_iterator 
         i=multiIndices.begin(); i != multiIndices.end(); i++)
    {
      if (i->order()==1)
        {
          subset->addDeriv(new CoordDeriv(i->firstOrderDirection()), 
                           VectorDeriv);
        }
      if (i->order()==0)
        {
          subset->addDeriv(MultipleDeriv(),
                           VectorDeriv);
        }
    }
  
  SUNDANCE_VERB_HIGH(tabs << "DiscreteFuncElement " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "DiscreteFuncElement " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));
  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

