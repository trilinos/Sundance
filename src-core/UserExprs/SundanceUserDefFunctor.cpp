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

#include "SundanceUserDefFunctor.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

UserDefFunctor::UserDefFunctor(const string& name)
  :  name_(name) 
{;}

void UserDefFunctor::getArgDerivIndices(const Array<int>& orders,
                                        Map<MultiSet<int>, int>& varArgDerivs,
                                        Map<MultiSet<int>, int>& constArgDerivs) const
{
  int n = numArgs();

  int k = 0;

  for (unsigned int o=0; o<orders.size(); o++)
    {
      int order = orders[o];

      if (order==0)
        {
          varArgDerivs.put(MultiSet<int>(), k++);
        }
      else if (order==1)
        {
          for (int i=0; i<n; i++)
            {
              varArgDerivs.put(makeMultiSet<int>(i), k++);
            }
        }
      else if (order==2)
        {
          for (int i=0; i<n; i++)
            {
              for (int j=0; j<n; j++)
                {
                  varArgDerivs.put(makeMultiSet<int>(i,j), k++);
                }
            }
        }
      else
        {
          TEST_FOR_EXCEPTION(order > 2 || order < 0, RuntimeError,
                             "order " << order << " not supported by functor " << name());
        }
    }
}



void UserDefFunctor::evalArgDerivs(int maxOrder,
                                   const Array<double>& args,
                                   Array<double>& argDerivs) const
{
  if (maxOrder==0)
    {
      argDerivs.resize(1);
      argDerivs[0] = eval0(args);
    }
  else if (maxOrder==1)
    {
      argDerivs.resize(args.size()+1);
      argDerivs[0] = eval1(args, &(argDerivs[1]));
    }
  else
    {
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "deriv order=" << maxOrder 
                         << " not supported for function " << name());
    }
}


double UserDefFunctor::eval0(const Array<double>& vars) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "eval0() function not supported for functor " << name());
  return 0.0; // -Wall
}

double UserDefFunctor::eval1(const Array<double>& vars, double* derivs) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError,
                     "eval1() function not supported for functor " << name());
  return 0.0; // -Wall
}
