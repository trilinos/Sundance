#include "SundanceStdMathFunctors.hpp"
#include "Teuchos_MPISession.hpp"

using namespace SundanceUtils;
using namespace Teuchos;

template <class F> bool functorTest(int nx, double tol)
{
  RefCountPtr<UnaryFunctor> f = rcp(new F());

  return f->test(nx, tol);
}


int main(int argc, void** argv)
{

  
  try
		{
      MPISession::init(&argc, &argv);

      int nx = 5;
      double tol = 1.0e-6;

      UnaryFunctor::checkResults() = true;

      UnaryFunctor::fdStep() = 1.0e-3;

      bool isOK = functorTest<StdReciprocal>(nx, tol);

      isOK = functorTest<StdExp>(nx, tol) && isOK ;

      isOK = functorTest<StdLog>(nx, tol) && isOK ;

      isOK = functorTest<StdSqrt>(nx, tol) && isOK ;

      /* trig functions */

      isOK = functorTest<StdSin>(nx, tol) && isOK ;

      isOK = functorTest<StdCos>(nx, tol) && isOK ;

      isOK = functorTest<StdTan>(nx, tol) && isOK ;

      isOK = functorTest<StdASin>(nx, tol) && isOK ;

      isOK = functorTest<StdACos>(nx, tol) && isOK ;

      isOK = functorTest<StdATan>(nx, tol) && isOK ;

      /* hyperbolic functions */

      isOK = functorTest<StdSinh>(nx, tol) && isOK ;

      isOK = functorTest<StdCosh>(nx, tol) && isOK ;

      isOK = functorTest<StdTanh>(nx, tol) && isOK ;

      isOK = functorTest<StdASinh>(nx, tol) && isOK ;

      isOK = functorTest<StdACosh>(nx, tol) && isOK ;

      isOK = functorTest<StdATanh>(nx, tol) && isOK ;

      cerr << "done with all tests!" << endl;

      if (isOK) 
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
