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

#ifndef SUNDANCE_COMBINATORIALUTILS_H
#define SUNDANCE_COMBINATORIALUTILS_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceMultiSet.hpp"

namespace SundanceUtils
{
  using namespace Teuchos;
  /**
   * Return partitions of an integer
   * @author Kevin Long
   */
  Array<Array<int> > partitionInteger(int n);

  /** 
   * Return compositions of an integer
   */
  Array<Array<Array<int> > > compositions(int n);

  /** */
  template <class T> inline
  Array<Array<Array<T> > > indexArrangements(const MultiSet<T>& mu,
                                             const Array<int>& k)
  {
    int nBins = k.size();
    
    int M = 0;
    for (int i=0; i<nBins; i++)
      {
        M += k[i];
      }
    
    Array<T> I;
    typename MultiSet<T>::const_iterator iter;
    for (iter=mu.begin(); iter!=mu.end(); iter++)
      {
        I.append(*iter);
      }

    Array<Array<Array<T> > > rtn;
    
    do
      {
        Array<Array<T> > bins(nBins);
        int count = 0;
        for (int i=0; i<nBins; i++)
          {
            for (int j=0; j<k[i]; j++)
              {
                bins[i].append(I[count++]);
              }
          }
        rtn.append(bins);
      }
    while (std::next_permutation(I.begin(), I.end()));
    return rtn;
  }

  /** 
   * Return the distinct arrangements of the given multiset into n bins
   */
  template <class T> inline
  Array<Array<Array<T> > > binnings(const MultiSet<T>& mu, int n)
  {
    int N = mu.size();
    Array<Array<int> > c = compositions(N)[n-1];
    Array<Array<Array<T> > > rtn;

    for (unsigned int i=0; i<c.size(); i++)
      {
        Array<Array<Array<T> > > a = indexArrangements(mu, c[i]);
        for (unsigned int j=0; j<a.size(); j++)
          {
            rtn.append(a[j]);
          }
      }
    return rtn;
  }


}

#endif  /* DOXYGEN_DEVELOPER_ONLY */   

#endif



