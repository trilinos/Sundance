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

#include "SundanceUnknownParameterElement.hpp"

#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

UnknownParameterElement
::UnknownParameterElement(const string& name,
                          const string& suffix,
                          int myIndex)
	: UnknownFuncElement(rcp(new UnknownFuncDataStub()), name, suffix, myIndex),
    SpatiallyConstantExpr()
{}

void UnknownParameterElement::findNonzeros(const EvalContext& context,
                                           const Set<MultiIndex>& multiIndices,
                                           const Set<MultiSet<int> >& inputActiveFuncIDs,
                                           bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for unknown parameter " 
                       << toString()
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

  if (evalPtIsZero())
    {
      SUNDANCE_VERB_MEDIUM(tabs << "UnknownParameterElement eval point is a zero expr");
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "UnknownParameterElement eval point is a nonzero expr");
    }

  /* Evaluate the function itself, i.e., the zeroth deriv of the function.
   * If this is a test function, or if we are doing a linear problem,
   * then we skip this step. */
  if (!regardFuncsAsConstant && !evalPtIsZero())
    {
      if (activeFuncIDs.contains(MultiSet<int>()))
        {
          Tabs tab1;
          SUNDANCE_VERB_MEDIUM(tab1 << "adding deriv {}");
          subset->addDeriv(MultipleDeriv(), ConstantDeriv);
          SUNDANCE_VERB_HIGH(tab1 << "sparsity subset is now "
                             << endl << *subset);
        }
      else
        {
          SUNDANCE_VERB_MEDIUM(tabs << "value of " << toString() << " not required");
        }
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "value of " << toString() << " not required");
    }
  
  /* If this function is one of the active variables, then
   * add the deriv wrt this func to the sparsity pattern */
  MultiSet<int> myFuncID;
  myFuncID.put(funcID());
  if (activeFuncIDs.contains(myFuncID))
    {
      subset->addDeriv(new FunctionalDeriv(this, MultiIndex()),
                       ConstantDeriv);
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "deriv wrt to " << toString() << " not required");
    }
  
  const Parameter* p
    = dynamic_cast<const Parameter*>(evalPt());
  if (p != 0)
    {
      Tabs tab1;
      SUNDANCE_VERB_MEDIUM(tabs << "UnknownParameterElement finding nonzero pattern for eval pt");
      p->findNonzeros(context, multiIndices, activeFuncIDs,
                      regardFuncsAsConstant);
    }
  

  SUNDANCE_VERB_HIGH(tabs << "unknown parameter " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "unknown parameter " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}


Set<MultipleDeriv> 
UnknownParameterElement::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) 
    {
      if (!evalPtIsZero()) rtn.put(MultipleDeriv());
    }
  else if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }

  return rtn;
}



Set<MultipleDeriv> 
UnknownParameterElement::internalFindV(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  return rtn;
}


Set<MultipleDeriv> 
UnknownParameterElement::internalFindC(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==1)
    {
      Deriv d = new FunctionalDeriv(this, MultiIndex());
      MultipleDeriv md;
      md.put(d);
      rtn.put(md);
    }
  return rtn.intersection(findR(order, context));
}



Evaluator* UnknownParameterElement
::createEvaluator(const EvaluatableExpr* expr,
                  const EvalContext& context) const 
{
  return SymbolicFuncElement::createEvaluator(expr, context);
}


const Parameter* UnknownParameterElement::parameterValue() const 
{
  const Parameter* p = dynamic_cast<const Parameter*>(evalPt());
  TEST_FOR_EXCEPTION(p==0, InternalError, 
                     "UnknownParameter evalPt() is not a Parameter");
  return p;
}

Parameter* UnknownParameterElement::parameterValue()  
{
  Parameter* p = dynamic_cast<Parameter*>(evalPt());
  TEST_FOR_EXCEPTION(p==0, InternalError, 
                     "UnknownParameter evalPt() is not a Parameter");
  return p;
}






XMLObject UnknownParameterElement::toXML() const 
{
	XMLObject rtn("UnknownParameterElement");
	rtn.addAttribute("name", name());
	return rtn;
}

