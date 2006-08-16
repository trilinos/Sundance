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

#include "SundanceUserDefOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUserDefOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;




UserDefOpEvaluator
::UserDefOpEvaluator(const UserDefOp* expr,
                     const EvalContext& context)
  : ChainRuleEvaluator(expr, context),
    argValueIndex_(expr->numChildren()),
    argValueIsConstant_(expr->numChildren()),
    functor_(expr->op()),
    maxOrder_(0),
    numVarArgDerivs_(0),
    numConstArgDerivs_(0)
{
  Tabs tab1;
  SUNDANCE_VERB_LOW(tab1 << "initializing user defined op evaluator for " 
                    << expr->toString());
  Array<int> orders = findRequiredOrders(expr, context);

  SUNDANCE_VERB_HIGH(tab1 << "setting arg deriv indices");
  
  
  /* Find the mapping from argument derivatives to indices in the 
   * functor's vector of return values */
  Map<MultiSet<int>, int> varArgDerivs;
  Map<MultiSet<int>, int> constArgDerivs;
  functor()->getArgDerivIndices(orders, varArgDerivs, constArgDerivs);
  numVarArgDerivs_ = varArgDerivs.size();
  numConstArgDerivs_ = constArgDerivs.size();
  typedef Map<MultiSet<int>, int>::const_iterator iter;
  for (iter i=varArgDerivs.begin(); i!=varArgDerivs.end(); i++)
    {
      Tabs tab2;
      SUNDANCE_VERB_EXTREME(tab2 << "variable arg deriv " << i->first << " will be at index "
                              << i->second);
      addVarArgDeriv(i->first, i->second);
    }
  
  for (iter i=constArgDerivs.begin(); i!=constArgDerivs.end(); i++)
    {
      Tabs tab2;
      SUNDANCE_VERB_EXTREME(tab2 << "constant arg deriv " << i->first << " will be at index "
                            << i->second);
      addConstArgDeriv(i->first, i->second);
    }

  /* Find the indices to the zeroth derivative of each argument */
  for (int i=0; i<expr->numChildren(); i++)
    {
      const SparsitySuperset* sArg = childSparsity(i);
      int numConst=0;
      int numVec=0;
      for (int j=0; j<sArg->numDerivs(); j++)
        {
          if (sArg->deriv(j).order() == 0) 
            {
              if (sArg->state(j)==VectorDeriv)
                {
                  argValueIndex_[i] = numVec;              
                }
              else
                {
                  argValueIndex_[i] = numConst;              
                }
              break;
            }
          if (sArg->state(j) == VectorDeriv) 
            {
              numVec++;
            }
          else
            {
              numConst++;
            }
        }
    }
  
  /* Call init() at the base class to set up chain rule evaluation */
  init(expr, context);
}


Array<int> UserDefOpEvaluator::findRequiredOrders(const ExprWithChildren* expr, 
                                                  const EvalContext& context)
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "finding required arg deriv orders");

  Set<int> orders;
  
  const Set<MultipleDeriv>& R = expr->findR(context);
  typedef Set<MultipleDeriv>::const_iterator iter;

  for (iter md=R.begin(); md!=R.end(); md++)
    {
      Tabs tab1;
      
      int N = md->order();
      if (N > maxOrder_) maxOrder_ = N;
      if (N==0) orders.put(N);
      for (int n=1; n<=N; n++)
        {
          const Set<MultiSet<int> >& QW = expr->findQ_W(n, context);
          for (Set<MultiSet<int> >::const_iterator q=QW.begin(); q!=QW.end(); q++)
            {
              orders.put(q->size());
            }
        }
    }
  SUNDANCE_VERB_HIGH(tab0 << "arg deriv orders=" << orders);
  return orders.elements();
}


void UserDefOpEvaluator
::evalArgDerivs(const EvalManager& mgr,
                const Array<RefCountPtr<Array<double> > >& constArgVals,
                const Array<RefCountPtr<Array<RefCountPtr<EvalVector> > > >& vArgVals,
                Array<double>& constArgDerivs,
                Array<RefCountPtr<EvalVector> >& varArgDerivs) const
{
  Array<double> argPoint(argValueIndex_.size());

  int numPoints = EvalManager::stack().vecSize();

  TEST_FOR_EXCEPTION(numPoints==0, InternalError,
                     "Empty vector detected in evalArgDerivs()");


  varArgDerivs.resize(numVarArgDerivs_);
  for (unsigned int i=0; i<varArgDerivs.size(); i++)
    {
      varArgDerivs[i] = mgr.popVector();
      varArgDerivs[i]->resize(numPoints);
      string str;
      if (i==0)
        {
          str = functor()->name() + "()";
        }
      else
        {
          str = functor()->name() + "_" + Teuchos::toString(i) + "()";
        }
      varArgDerivs[i]->setString(str);
    }
  
  Array<double> argDerivsAtPoint(varArgDerivs.size());

  for (int p=0; p<numPoints; p++)
    {
      Tabs tab1;
      /* get the argument value at the p-th point */
      for (unsigned int q=0; q<argValueIndex_.size(); q++)
        {
          Tabs tab2;
          if (argValueIsConstant_[q]) 
            {
              argPoint[q] = (*(constArgVals[q]))[argValueIndex_[q]];
            }
          else
            {
              argPoint[q] = (*(vArgVals[q]))[argValueIndex_[q]]->start()[p];
            }
        }

      /* evaluate the function at the point */      
      functor()->evalArgDerivs(maxOrder_, argPoint, argDerivsAtPoint);


      /* Put the results of the functor call into the arg deriv array */
      for (unsigned int i=0; i<varArgDerivs.size(); i++)
        {
          varArgDerivs[i]->start()[p] = argDerivsAtPoint[i];
        }
    }
}

