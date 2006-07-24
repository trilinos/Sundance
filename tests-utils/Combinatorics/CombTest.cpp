#include "SundanceCombinatorialUtils.hpp"
#include "Teuchos_MPISession.hpp"


using namespace SundanceUtils;
using namespace Teuchos;




int main(int argc, void** argv)
{
  try
		{
      MPISession::init(&argc, &argv);


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

          /*

          for (unsigned int i=0; i<x.size(); i++)
            {
              for (unsigned int j=0; j<x[i].size(); j++)
                {
                  int sum=0;
                  for (unsigned int k=0; k<x[i][j].size(); k++)
                    {
                      sum += x[i][j][k];
                    }
                  if (sum != n) bad = true;
                  cout << "arrangements(" << x[i][j] << " = " 
                       << indexArrangements(x[i][j]) << endl;
                }
            }

          */
        }


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

  MPISession::finalize();
}
