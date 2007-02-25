#include "SundanceCombinatorialUtils.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace SundanceUtils;
using namespace Teuchos;



#define TEST_MS(x) \
  {\
    Set<MultiSet<int> > subs = multisetSubsets(x);\
    write(x, subs);                                                \
    Set<MultiSet<MultiSet<int> > > parts = multisetPartitions(x);\
    write(x, parts);                                                 \
    Array<Array<MultiSet<int> > > comps = multisetCompositions(x);\
    write(x, comps);                                               \
  }



void write(const MultiSet<int>& x, 
           const Set<MultiSet<MultiSet<int> > >& y)
{
  cout << "---- Partitions of " << x << " ------------------"
       << endl;
  for (Set<MultiSet<MultiSet<int> > >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      cout << *i << endl;
    }
}

void write(const MultiSet<int>& x, 
           const Array<Array<MultiSet<int> > >& y)
{
  cout << "---- Compositions of " << x << " ------------------"
       << endl;
  for (unsigned int i=0; i<y.size(); i++)
    {
      cout << y[i] << endl;
    }
}

void write(const MultiSet<int>& x, 
           const Set<MultiSet<int> >& y)
{
  cout << "---- Subsets of " << x << " ------------------"
       << endl;
  for (Set<MultiSet<int> >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      cout << *i << endl;
    }
}



int main(int argc, char** argv)
{
  try
		{
      GlobalMPISession session(&argc, &argv);


      bool bad = false;

      for (int n=1; n<=4; n++)
        {
          Array<Array<Array<int> > > c = compositions(n);
          cout << "N=" << n << " compositions=" << c << endl;

          MultiSet<int> mu;
          for (int m=1; m<=n; m++)
            {
              mu.put(m);

            }
          for (int m=1; m<=n; m++)
            {
              Array<Array<Array<int> > > b = binnings(mu, m);
              cout << "binnings = " << b << endl;
            }

          cout << "--------- non-neg compositions" << endl;
          for (int m=1; m<=n; m++)
            {
              for (int k=1; k<=n; k++)
                {
                  Array<Array<int> > a = nonNegCompositions(m, k);
                  cout << m << " " << k << " " << endl;
                  for (unsigned int l=0; l<a.size(); l++)
                    {
                      cout << "         " << a[l] << endl;
                    }
                }
            }
          
          cout << "-------- index combs ---- " << endl;
          Array<int> s = tuple(2,3,2);
          Array<Array<int> > C = indexCombinations(s);
          for (unsigned int m=0; m<C.size(); m++)
            {
              cout << C[m] << endl;
            }
        }

      cout << "--------- index tuples ----------------" << endl;

      Array<Array<int> > x = distinctIndexTuples(2, 6);

      cout << "num choices = " << x.size() << endl;

      for (unsigned int i=0; i<x.size(); i++) 
        {
          if ((i % 5)==0) cout << endl;
          cout << x[i] << endl;
        }
      
#ifdef BLAH
      TEST_MS(makeMultiSet(1));
      TEST_MS(makeMultiSet(1, 1));
      TEST_MS(makeMultiSet(1, 2));
      TEST_MS(makeMultiSet(1, 1, 2));
      TEST_MS(makeMultiSet(1, 1, 2, 2));
      TEST_MS(makeMultiSet(1, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3, 3));
#endif BLAH

      if (!bad) 
        {
          cerr << "all tests PASSED" << endl;
        }
      else
        {
          cerr << "a test has FAILED" << endl;
        }
    }
	catch(std::exception& e)
		{
      cerr << "detected exception " << e.what() << endl;
		}

  
}
