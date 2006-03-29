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

#include "SundanceMultipleDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceFunctionalDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

MultipleDeriv::MultipleDeriv()
  : MultiSet<Deriv>()
{}

int MultipleDeriv::spatialOrder() const 
{
  int rtn = 0;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
    {
      if (i->isCoordDeriv())
        {
          rtn += 1;
        }
    }
  return rtn;
}

MultiIndex MultipleDeriv::spatialDeriv() const
{
  MultiIndex rtn;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
    {
      if (i->isCoordDeriv())
        {
          int d = i->coordDeriv()->dir();
          rtn[d] += 1;
        }
    }
  return rtn;
}

MultiSet<int> MultipleDeriv::funcIDs() const
{
  MultiSet<int> rtn;
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
    {
      if (i->isFunctionalDeriv())
        {
          int f = i->funcDeriv()->funcID();
          rtn.put(f);
        }
      TEST_FOR_EXCEPTION(!i->isFunctionalDeriv(), RuntimeError,
                         "MultipleDeriv::funcIDs() found spatial deriv");
    }
  return rtn;
}

MultipleDeriv MultipleDeriv::product(const MultipleDeriv& other) const 
{
  MultipleDeriv rtn;
  
  for (MultipleDeriv::const_iterator i=this->begin(); i!=this->end(); i++)
    {
      rtn.put(*i);
    }
  for (MultipleDeriv::const_iterator i=other.begin(); i!=other.end(); i++)
    {
      rtn.put(*i);
    }
  return rtn;
}

MultipleDeriv MultipleDeriv::factorOutDeriv(const Deriv& x) const
{
  MultipleDeriv rtn;
  MultipleDeriv::const_iterator i;
  bool gotIt = false;
  for (i=this->begin(); i!=this->end(); i++)
    {
      if (!gotIt && *i==x) {gotIt=true; continue;}
      rtn.put(*i);
    }
  return rtn;
}



void MultipleDeriv
::productRulePermutations(ProductRulePerms& perms) const 
{
  int N = order();

  if (N==0)
    {
      MultipleDeriv md0;
      DerivPair p(md0, md0);
      perms.put(p, 1);
      return;
    }

  int p2 = pow2(N);

  for (int i=0; i<p2; i++)
    {
      MultipleDeriv left;
      MultipleDeriv right;
      Array<int> bits = bitsOfAnInteger(i, N);
      int j=0; 
      MultipleDeriv::const_iterator iter;
      for (iter=this->begin(); iter != this->end(); iter++, j++)
        {
          if (bits[j]==true)
            {
              left.put(*iter);
            }
          else
            {
              right.put(*iter);
            }
        }
      DerivPair p(left, right);
      if (!perms.containsKey(p))
        {
          perms.put(p, 1);
        }
      else
        {
          int count = perms.get(p);
          perms.put(p, count+1);
        }
    }
}

Array<int> MultipleDeriv::bitsOfAnInteger(int x, int n)
{
  TEST_FOR_EXCEPTION(x >= pow2(n), InternalError,
                     "Invalid input to MultipleDeriv::bitsOfX");
                     
  Array<int> rtn(n);

  int r = x;
  for (int b=n-1; b>=0; b--)
    {
      rtn[b] = r/pow2(b);
      r = r - rtn[b]*pow2(b);
    }
  return rtn;
}

int MultipleDeriv::pow2(int n)
{
  static Array<int> p2(1,1);

  if ((unsigned int) n >= p2.size())
    {
      int oldN = p2.size(); 
      for (int i=oldN; i<=n; i++) p2.push_back(p2[i-1]*2);
    }
  
  return p2[n];
}

