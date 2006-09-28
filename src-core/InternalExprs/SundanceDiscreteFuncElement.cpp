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
	: EvaluatableExpr(), 
    FuncElementBase(name, suffix),
    commonData_(data),
    myIndex_(myIndex),
    miSet_()
{}


RefCountPtr<Array<Set<MultipleDeriv> > > DiscreteFuncElement
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "DFE::internalDetermineR() for "
                     << toString());
  Array<Set<MultipleDeriv> > RIn = RInput;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
    {
      const MultiIndex& mi = *i;
      int order = mi.order();
      if (order==0) RIn[0].put(MultipleDeriv());
      if (order==1) RIn[1].put(MultipleDeriv(new CoordDeriv(mi.firstOrderDirection())));
    }

  return EvaluatableExpr::internalDetermineR(context, RIn);
}


Set<MultipleDeriv> 
DiscreteFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
    {
      if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
    }
  if (order==1)
    {
      for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
        {
          const MultiIndex& mi = *i;
          int diffOrder = mi.order();
          if (diffOrder==1) 
            rtn.put(MultipleDeriv(new CoordDeriv(mi.firstOrderDirection())));
        }
    }

  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "DFE::internalFindV(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
    {
      if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
    }
  if (order==1)
    {
      for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
        {
          const MultiIndex& mi = *i;
          int diffOrder = mi.order();
          if (diffOrder==1) 
            rtn.put(MultipleDeriv(new CoordDeriv(mi.firstOrderDirection())));
        }
    }
  
  rtn = rtn.intersection(findR(order, context));
  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  return rtn;
}

void DiscreteFuncElement::addMultiIndex(const MultiIndex& newMi) const
{
  miSet_.put(newMi);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}



