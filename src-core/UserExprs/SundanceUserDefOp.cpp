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
#include "SundanceUserDefOpElement.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefFunctorElement.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore;
using namespace Teuchos;
using namespace TSFExtended;

UserDefOp::UserDefOp(const Expr& args,
                     const RefCountPtr<const UserDefFunctor>& functor)
  : ListExpr()
{
  int nElems = functor->rangeDim();
  Array<RefCountPtr<ScalarExpr> > scalarArgs = getScalarArgs(args);

  RefCountPtr<SundanceUtils::Map<EvalContext, RefCountPtr<const UserDefOpCommonEvaluator> > > 
    commonEvaluatorsMap = rcp(new SundanceUtils::Map<EvalContext, RefCountPtr<const UserDefOpCommonEvaluator> >());
  
  for (int i=0; i<nElems; i++)
    {
      RefCountPtr<UserDefFunctorElement> e 
        = rcp(new UserDefFunctorElement(functor, i));
      Expr elem = new UserDefOpElement(scalarArgs, commonEvaluatorsMap, e);
      append(elem);
    }
}



Array<RefCountPtr<ScalarExpr> > UserDefOp::getScalarArgs(const Expr& args)
{
  Expr fargs = args.flatten();
  Array<RefCountPtr<ScalarExpr> > sargs(fargs.size());
  
  for (unsigned int i=0; i<fargs.size(); i++)
    {
      sargs[i] = rcp_dynamic_cast<ScalarExpr>(fargs[i].ptr());
    }
  return sargs;
}

