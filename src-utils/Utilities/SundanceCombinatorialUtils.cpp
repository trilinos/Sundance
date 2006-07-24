#include "SundanceCombinatorialUtils.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceTabs.hpp"
#include <algorithm>
#include <iterator>
#include <iostream>


namespace SundanceUtils
{
  using namespace Teuchos;

  Array<Array<int> > partitionInteger(int n)
  {
    static Array<Array<Array<int> > > rtn 
      = tuple(
              tuple(tuple(1)), 
              tuple(tuple(2), tuple(1, 1)),
              tuple(tuple(3), tuple(2, 1), tuple(1, 1, 1)),
              tuple(tuple(4), tuple(3, 1), tuple(2, 2), tuple(2, 1, 1), tuple(1,1,1,1)));
    TEST_FOR_EXCEPTION(n<1 || n>4, RuntimeError, 
                       "case n=" << n << " not implemented in partitionInteger()");
    return rtn[n-1];
  }






  Array<Array<Array<int> > > compositions(int n)
  {
    Array<Array<Array<int> > > q(n);

    Array<Array<int> > x = partitionInteger(n);

    Array<Array<int> > p;
    for (unsigned int m=0; m<x.size(); m++)
      {
        Array<int> tmp;
        vector<int> y = x[m];
        std::sort(y.begin(), y.end());
        tmp.resize(y.size());
        copy(y.begin(), y.end(), tmp.begin());
        p.append(tmp);
        while (std::next_permutation(y.begin(), y.end())) 
          { 
            tmp.resize(y.size());
            copy(y.begin(), y.end(), tmp.begin());
            p.append(tmp);
          }
      }

    for (unsigned int i=0; i<p.size(); i++)
      {
        q[p[i].size()-1].append(p[i]);
      }

    return q;
  }



  Array<Array<Array<int> > > binnings(const MultiSet<int>& mu, int n)
  {
    int N = mu.size();
    Array<Array<int> > c = compositions(N)[n-1];
    Array<Array<Array<int> > > rtn;

    for (unsigned int i=0; i<c.size(); i++)
      {
        Array<Array<Array<int> > > a = indexArrangements(mu, c[i]);
        for (unsigned int j=0; j<a.size(); j++)
          {
            rtn.append(a[j]);
          }
      }
    return rtn;
  }



  Array<Array<Array<int> > > indexArrangements(const MultiSet<int>& mu,
                                               const Array<int> & k) 
  {
    int nBins = k.size();
    
    int M = 0;
    for (int i=0; i<nBins; i++)
      {
        M += k[i];
      }

    Array<int> I;
    for (MultiSet<int>::const_iterator iter=mu.begin(); iter!=mu.end(); iter++)
      {
        I.append(*iter);
      }

    Array<Array<Array<int> > > rtn;
    
    do
      {
        Array<Array<int> > bins(nBins);
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
}



	





