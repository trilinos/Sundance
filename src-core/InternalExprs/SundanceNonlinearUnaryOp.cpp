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

#include "SundanceNonlinearUnaryOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

NonlinearUnaryOp::NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg,
                                   const RefCountPtr<UnaryFunctor>& op)
  : UnaryExpr(arg), op_(op)
{
  if (isEvaluatable(arg.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          if (evaluatableArg()->orderOfSpatialDependency(d) != 0) 
            {
              setOrderOfDependency(d, -1);
            }
          else
            {
              setOrderOfDependency(d, 0);
            }
        }

      const Set<MultiSet<int> >& argFuncs = evaluatableArg()->funcIDSet();
      typedef Set<MultiSet<int> >::const_iterator iter;

      for (iter i=argFuncs.begin(); i != argFuncs.end(); i++)
        {
          const MultiSet<int>& f1 = *i;
          addFuncIDCombo(f1);
          for (iter j=argFuncs.begin(); j != argFuncs.end(); j++)
            {
              const MultiSet<int>& f2 = *j;
              MultiSet<int> f12 = f1.merge(f2);
              if (f1.size()+f2.size() > maxFuncDiffOrder()) continue;
              addFuncIDCombo(f12);
              for (iter k=argFuncs.begin(); k != argFuncs.end(); k++)
                {
                  const MultiSet<int>& f3 = *k;
                  
                  if (f1.size()+f2.size()+f3.size() > maxFuncDiffOrder()) 
                    continue;
                  addFuncIDCombo(f12.merge(f3));
                }
            }
        }
    }
}

Set<MultiIndex> NonlinearUnaryOp
::argMultiIndices(const Set<MultiIndex>& miSet) const 
{
  Set<MultiIndex> rtn = miSet;
  rtn.put(MultiIndex());
  return rtn;
}


Set<MultiSet<int> > 
NonlinearUnaryOp::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                                 int inputMaxOrder) const 
{
  typedef Set<MultiSet<int> >::const_iterator iter;

  //  const Set<int>& argFuncs = evaluatableArg()->funcDependencies();
  const Set<MultiSet<int> >& argFuncIDs = evaluatableArg()->funcIDSet();

  Set<MultiSet<int> > rtn;
  if (activeFuncIDs.contains(MultiSet<int>())) rtn.put(MultiSet<int>());

  int maxOrder = 0;
  if (true/*inputMaxOrder < 0*/)
    {
      for (iter i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
        {
          const MultiSet<int>& f1 = *i;
          maxOrder = max(maxOrder, (int) f1.size());
        }
      maxOrder = min(maxOrder, (int) maxFuncDiffOrder());
    }
  else
    {
      maxOrder = inputMaxOrder;
    }

  for (iter i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& f1 = *i;
      if (((int)f1.size()) > maxOrder) continue;
      if (!argFuncIDs.contains(f1)) continue;
      rtn.put(f1);
      for (iter j=activeFuncIDs.begin(); j != activeFuncIDs.end(); j++)
        {
          const MultiSet<int>& f2 = *j;
          if ( ((int) (f1.size()+f2.size())) > maxOrder) continue;
          MultiSet<int> f12 = f1.merge(f2);
          if (!argFuncIDs.contains(f12)) continue;
          rtn.put(f12);
          for (iter k=activeFuncIDs.begin(); k != activeFuncIDs.end(); k++)
            {
              const MultiSet<int>& f3 = *k;
              if ( ((int) (f1.size()+f2.size() +f3.size())) > maxOrder) continue;
              MultiSet<int> f123 = f12.merge(f3);
              if (!argFuncIDs.contains(f123)) continue;
              rtn.put(f123);
            }
        }
    }
  /* if we're going to compute anything, we need to compute the argument */
  if (rtn.size() > 0) rtn.put(MultiSet<int>());


  return rtn;
}

void NonlinearUnaryOp::findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    const Set<MultiSet<int> >& inputActiveFuncIDs,
                                    bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for unary nonlinear op " 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  Set<MultiSet<int> > activeFuncIDs = filterActiveFuncs(inputActiveFuncIDs);


  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  addActiveFuncs(context, activeFuncIDs);
  RefCountPtr<SparsitySubset> subset 
    = sparsitySubset(context, multiIndices, activeFuncIDs, false);

  Set<MultiSet<int> > childFuncIDs = argActiveFuncs(activeFuncIDs, -1);


  
  SUNDANCE_VERB_MEDIUM(tabs << "NonlinUnaryOp arg funcID set is " << childFuncIDs);

  int maxMiOrder = maxOrder(multiIndices);
  int maxDiffOrder = context.topLevelDiffOrder() + maxMiOrder;

  Set<MultiIndex> argMI = argMultiIndices(multiIndices);

  evaluatableArg()->findNonzeros(context, argMI,
                                 childFuncIDs,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsitySubset 
    = evaluatableArg()->sparsitySubset(context, argMI, childFuncIDs, true);

  SUNDANCE_VERB_MEDIUM(tabs << "NonlinUnaryOp arg sparsity subset is " 
                       << endl << *argSparsitySubset);

  SUNDANCE_VERB_MEDIUM(tabs << "active func IDs are " 
                       << endl << activeFuncIDs);
  

  /* Determine whether the argument is a constant or a vector */
  DerivState argState = VectorDeriv;
  for (int i=0; i<argSparsitySubset->numDerivs(); i++)
    {
      if (argSparsitySubset->deriv(i).order()==0)
        {
          argState = argSparsitySubset->state(i);
          break;
        }
    }

  for (int i=0; i<argSparsitySubset->numDerivs(); i++)
    {
      const MultipleDeriv& d = argSparsitySubset->deriv(i);
      if (!activeFuncIDs.contains(d.funcIDs())) 
        {
          Tabs tab3;
          SUNDANCE_VERB_MEDIUM(tab3 << "skipping inactive deriv" << d);
          continue;
        }
      subset->addDeriv(d, argState);
    }

  for (int i=0; i<argSparsitySubset->numDerivs(); i++)
    {
      for (int j=0; j<argSparsitySubset->numDerivs(); j++)
        {
          MultipleDeriv product 
            = argSparsitySubset->deriv(i).product(argSparsitySubset->deriv(j));
          if (product.order() > maxDiffOrder || product.spatialOrder() > maxMiOrder) 
            {
              Tabs tab3;
              SUNDANCE_VERB_MEDIUM(tab3 << "skipping deriv" << product
                                   << ": order too high");
            }
          if (!activeFuncIDs.contains(product.funcIDs())) 
            {
              Tabs tab3;
              SUNDANCE_VERB_MEDIUM(tab3 << "skipping inactive deriv" << product);
              continue;
            }
          subset->addDeriv(product, argState);
        }
    }



  SUNDANCE_VERB_HIGH(tabs << "nonlinear op  " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "nonlinear op " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}


ostream& NonlinearUnaryOp::toText(ostream& os, bool paren) const 
{
  os << op_->name() << "(" << arg().toString() << ")";
  return os;
}

ostream& NonlinearUnaryOp::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject NonlinearUnaryOp::toXML() const
{
  XMLObject rtn("NonlinearUnaryOp");
  rtn.addChild(arg().toXML());
  rtn.addAttribute("op", op_->name());
  return rtn;
}

