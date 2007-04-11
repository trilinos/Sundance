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

#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


SymbolicFuncElement::SymbolicFuncElement(const string& name,
  const string& suffix,
  int commonFuncID,
  int myIndex)
	: EvaluatableExpr(), FuncElementBase(name, suffix, commonFuncID),
    evalPt_(),
    evalPtDerivSetIndices_(),
    myIndex_(myIndex)
{}

void SymbolicFuncElement::registerSpatialDerivs(const EvalContext& context, 
                                                const Set<MultiIndex>& miSet) const
{
  evalPt()->registerSpatialDerivs(context, miSet);
  EvaluatableExpr::registerSpatialDerivs(context, miSet);
}

Set<MultipleDeriv> 
SymbolicFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "SFE::internalFindW(order=" << order << ") for "
                     << toString());

  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_VERB_HIGH(tab1 << "findW() for eval point");
    evalPt()->findW(order, context);
  }

  if (order==0) 
    {
      Tabs tab1;
      if (!evalPtIsZero()) 
        {
          SUNDANCE_VERB_EXTREME(tab1 << "value of " << toString() << " is nonzero" );
          rtn.put(MultipleDeriv());
        }
      else
        {
          SUNDANCE_VERB_EXTREME( tab1 << "value of " << toString() << " is zero" );
        }
    }
  else if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }
  
  SUNDANCE_VERB_HIGH( tab << "SFE: W[" << order << "] = " << rtn );

  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "SFE::internalFindV(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_VERB_HIGH(tab1 << "findV() for eval point");
    evalPt()->findV(order, context);
  }

  if (order==0) 
    {
      if (!evalPtIsZero()) rtn.put(MultipleDeriv());
    }

  SUNDANCE_VERB_EXTREME( tab << "SFE: V = " << rtn );
  SUNDANCE_VERB_EXTREME( tab << "SFE: R = " << findR(order, context) );
  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_VERB_HIGH( tab << "SFE: V[" << order << "] = " << rtn );
  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "SFE::internalFindV(order=" << order << ") for "
                     << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_VERB_HIGH(tab1 << "findC() for eval point");
    evalPt()->findC(order, context);
  }

  if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }

  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_VERB_HIGH( tab << "SFE: C[" << order << "] = " << rtn );
  return rtn;
}


RefCountPtr<Array<Set<MultipleDeriv> > > SymbolicFuncElement
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "SFE::internalDetermineR() for "
                     << toString());
  {
    Tabs tab1;
    SUNDANCE_VERB_HIGH(tab1 << "determineR() for eval point");
    evalPt()->determineR(context, RInput);
  }

  return EvaluatableExpr::internalDetermineR(context, RInput);
}

bool SymbolicFuncElement::evalPtIsZero() const
{
  TEST_FOR_EXCEPTION(evalPt_.get()==0, InternalError,
                     "null evaluation point in SymbolicFuncElement::evalPtIsZero()");
  bool isZero = 0 != dynamic_cast<const ZeroExpr*>(evalPt());
  bool isTest = 0 != dynamic_cast<const TestFuncElement*>(this);
  return isZero || isTest;
}

void SymbolicFuncElement::substituteZero() const 
{
  evalPt_ = rcp(new ZeroExpr());
}

void SymbolicFuncElement
::substituteFunction(const RefCountPtr<DiscreteFuncElement>& u0) const
{
  evalPt_ = u0;
}





