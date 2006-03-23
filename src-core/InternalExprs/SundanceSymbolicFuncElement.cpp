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
                                         int myIndex)
	: EvaluatableExpr(), FuncElementBase(name, suffix),
    evalPt_(),
    evalPtDerivSetIndices_(),
    myIndex_(myIndex)
{
  /* I have nonzero functional deriv of zero order, and of first order wrt 
   * myself */
  int fid = funcID();
  MultiSet<int> derivWrtMe;
  derivWrtMe.put(fid);
  addFuncIDCombo(derivWrtMe);
  addFuncIDCombo(MultiSet<int>());
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

void SymbolicFuncElement
::findNonzeros(const EvalContext& context,
               const Set<MultiIndex>& multiIndices,
               const Set<MultiSet<int> >& inputActiveFuncIDs,
               bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for symbolic func " 
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
      SUNDANCE_VERB_MEDIUM(tabs << "SymbolicFuncElement eval point is a zero expr");
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "SymbolicFuncElement eval point is a nonzero expr");
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
          subset->addDeriv(MultipleDeriv(), VectorDeriv);
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
  
  const DiscreteFuncElement* df 
    = dynamic_cast<const DiscreteFuncElement*>(evalPt());
  if (df != 0)
    {
      Tabs tab1;
      SUNDANCE_VERB_MEDIUM(tabs << "SymbolicFuncElement finding nonzero pattern for eval pt");
      df->findNonzeros(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant);
    }
  

  SUNDANCE_VERB_HIGH(tabs << "symbolic func " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "symbolic func " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}






