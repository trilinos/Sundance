//@HEADER@

//@HEADER@


#include <cstdlib>
#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearCombinationTester.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif



using namespace Teuchos;
using namespace std;
using namespace Playa;
using namespace PlayaOps;


int main(int argc, char *argv[]) 
{
  int stat = 0;

  try
    {
      GlobalMPISession session(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int nLocalRows = 10;
      
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;

      LinearCombinationTester<double> tester(nLocalRows,
                                             onProcDensity,
                                             offProcDensity,
                                             type,
                                             TestSpecifier<double>(true, 1.0e-13, 1.0e-10));

      
      bool allPass = tester.runAllTests();
      if (!allPass) stat = -1;
    }
  catch(std::exception& e)
    {
      stat = -1;
      std::cerr << "Caught exception: " << e.what() << std::endl;
    }
  return stat;
}



