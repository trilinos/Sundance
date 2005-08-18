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

#include "SundanceUserDefOp.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

UserDefOp::UserDefOp(const Expr& arg,
                     const RefCountPtr<UserDefFunctor>& op)
  : ExprWithChildren(getScalarArgs(arg)), op_(op)
{
  Set<MultiSet<int> > funcCombos;

  for (int i=0; i<numChildren(); i++)
    {
      if (isEvaluatable(evaluatableChild(i)))
        {
          for (int d=0; d<MultiIndex::maxDim(); d++) 
            {
              if (evaluatableChild(i)->orderOfSpatialDependency(d) != 0) 
                {
                  setOrderOfDependency(d, -1);
                }
              else
                {
                  setOrderOfDependency(d, 0);
                }
            }
          funcCombos.merge(evaluatableChild(i)->funcIDSet());
        }
    }
  
  typedef Set<MultiSet<int> >::const_iterator iter;

  for (iter i=funcCombos.begin(); i != funcCombos.end(); i++)
    {
      const MultiSet<int>& f1 = *i;
      for (iter j=funcCombos.begin(); j != funcCombos.end(); j++)
        {
          const MultiSet<int>& f2 = *j;
          MultiSet<int> f12 = f1.merge(f2);
          for (iter k=funcCombos.begin(); k != funcCombos.end(); k++)
            {
              const MultiSet<int>& f3 = *k;
              
              if (f1.size()+f2.size()+f3.size() > maxFuncDiffOrder()) 
                continue;
              addFuncIDCombo(f12.merge(f3));
            }
        }
    }

  
}

Set<MultiSet<int> > 
UserDefOp::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs) const 
{
  typedef Set<MultiSet<int> >::const_iterator iter;

  Set<MultiSet<int> > rtn;
  for (iter i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& f1 = *i;
      for (iter j=activeFuncIDs.begin(); j != activeFuncIDs.end(); j++)
        {
          const MultiSet<int>& f2 = *j;
          MultiSet<int> f12 = f1.merge(f2);
          for (iter k=activeFuncIDs.begin(); k != activeFuncIDs.end(); k++)
            {
              const MultiSet<int>& f3 = *k;
              
              if (f1.size()+f2.size()+f3.size() > maxFuncDiffOrder()) 
                continue;
              rtn.put(f12.merge(f3));
            }
        }
    }
  return rtn;
}


ostream& UserDefOp::toText(ostream& os, bool paren) const 
{
  os << op_->name() << "(";
  for (int i=0; i<numChildren(); i++)
    {
      os << child(i).toString();
      if (i < numChildren()-1) os << ",";
    }
  os << ")";
  return os;
}

ostream& UserDefOp::toLatex(ostream& os, bool paren) const 
{
  return toText(os, paren);
}

XMLObject UserDefOp::toXML() const
{
  XMLObject rtn("UserDefOp");
  XMLObject args("Arguments");
  for (int i=0; i<numChildren(); i++)
    {
      args.addChild(child(i).toXML());
    }
  rtn.addChild(args);
  rtn.addAttribute("op", op_->name());
  return rtn;
}

void UserDefOp::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                             const Set<MultiSet<int> >& activeFuncIDs,
                             bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for user-defined "
                       "nonlinear op" << toString() 
                       << " subject to multiindices " << multiIndices); 
  SUNDANCE_VERB_MEDIUM(tabs << "active funcs are " << activeFuncIDs);



  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  Set<MultiSet<int> > childFuncIDs = argActiveFuncs(activeFuncIDs);
  SUNDANCE_VERB_HIGH(tabs << "active funcs for arg are: " << childFuncIDs);

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs);

  for (int i=0; i<numChildren(); i++)
    {
      evaluatableChild(i)->findNonzeros(context, multiIndices,
                                        childFuncIDs,
                                        regardFuncsAsConstant);

      RefCountPtr<SparsitySubset> argSparsitySubset 
        = evaluatableChild(i)->sparsitySubset(context, multiIndices, childFuncIDs);

      for (int j=0; j<argSparsitySubset->numDerivs(); j++)
        {
          subset->addDeriv(argSparsitySubset->deriv(j), VectorDeriv);
        }
    }

  SUNDANCE_VERB_HIGH(tabs << "user-defined op " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "user-defined op " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}


Array<RefCountPtr<ScalarExpr> > UserDefOp::getScalarArgs(const Expr& args)
{
  Expr fargs = args.flatten();
  Array<RefCountPtr<ScalarExpr> > sargs(fargs.size());
  
  for (unsigned int i=0; i<fargs.size(); i++)
    {
      sargs[i] = rcp_dynamic_cast<ScalarExpr>(args[i].ptr());
    }
  return sargs;
}
