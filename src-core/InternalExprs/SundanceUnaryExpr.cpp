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

#include "SundanceUnaryExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


UnaryExpr::UnaryExpr(const RefCountPtr<ScalarExpr>& arg)
	: ExprWithChildren(tuple(arg)), allActiveFuncs_()
{}


void UnaryExpr::addActiveFuncs(const EvalContext& context,
                               const Set<MultiSet<int> >& activeFuncIDs) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "UnaryExpr " << toString() << " adding active funcs "
                     << activeFuncIDs);


  typedef Set<MultiSet<int> >::const_iterator iter;
  Set<MultiSet<int> > funcs = evaluatableArg()->funcIDSet();
  Set<MultiSet<int> > myActiveFuncs = activeFuncIDs.intersection(funcs);  

  if (allActiveFuncs_.containsKey(context))
    {
      allActiveFuncs_[context].merge(myActiveFuncs);        
    }
  else
    {
      allActiveFuncs_.put(context, myActiveFuncs);
    }
  
  SUNDANCE_VERB_HIGH(tabs << "UnaryExpr " << toString() << " added active funcs "
                     << myActiveFuncs);

}

const Set<MultiSet<int> >& UnaryExpr::getActiveFuncs(const EvalContext& context) const 
{
  TEST_FOR_EXCEPTION(!allActiveFuncs_.containsKey(context), 
                     InternalError,
                     "context " << context << " does not exist in UnaryExpr::getActiveFuncs()");
  return allActiveFuncs_.get(context);
}

Set<MultiSet<int> > 
UnaryExpr::argActiveFuncs(const Set<MultiSet<int> >& activeFuncID, 
                          int maxOrder) const
{
  typedef Set<MultiSet<int> >::const_iterator iter;
  Set<MultiSet<int> > rtn;

  Set<MultiSet<int> > funcs = evaluatableArg()->funcIDSet();


  for (iter i=activeFuncID.begin(); i!=activeFuncID.end(); i++)
    {
      if (maxOrder < 0 || (maxOrder>=0 && ((int) i->size()) <= maxOrder)) 
        {
          if (funcs.contains(*i)) rtn.put(*i);
        }
    }
  return rtn;
}
