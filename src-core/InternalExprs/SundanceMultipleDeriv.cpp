/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMultipleDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

MultipleDeriv::MultipleDeriv()
  : MultiSet<Deriv>()
{}

void MultipleDeriv::productRulePermutations(Array<MultipleDeriv>& left,
                                            Array<MultipleDeriv>& right) const
{
  int N = order();

  if (N==0)
    {
      left = Array<MultipleDeriv>(1);
      right = Array<MultipleDeriv>(1);
      MultipleDeriv d0;
      left[0] = d0;
      right[0] = d0;
      return;
    }

  int p2 = pow2(N);

  left = Array<MultipleDeriv>(p2);
  right = Array<MultipleDeriv>(p2);

  std::multiset<Deriv>::const_iterator iter;

  for (int i=0; i<p2; i++)
    {
      Array<int> bits = bitsOfAnInteger(i, N);
      int j=0; 
      for (iter=begin(); iter != end(); iter++, j++)
        {
          if (bits[j]==true)
            {
              left[i].put(*iter);
            }
          else
            {
              right[i].put(*iter);
            }
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

  if (n >= p2.size())
    {
      int oldN = p2.size(); 
      for (int i=oldN; i<=n; i++) p2.push_back(p2[i-1]*2);
    }
  
  return p2[n];
}

