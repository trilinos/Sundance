#include "Sundance.hpp"

#define DO_TEST(testName)\
  { \
    bool testName(); \
    bool pass = testName(); \
    if (pass) \
    {\
      numPass++; \
      Out::root() << "test " #testName " PASSED!" << endl; \
    } \
    else \
    { \
      numFail++; \
      Out::root() << "test " #testName " FAILED!" << endl; \
    }\
  }


int main(int argc, char** argv)
{
  bool allPass = true;

  try
  {
    Sundance::init(&argc, &argv);

    int numPass = 0;
    int numFail = 0;

    DO_TEST(BlockStochPoissonTest1D);
    DO_TEST(HighOrderPoisson2D);
    DO_TEST(HighOrderPoissonBernstein2D);
    DO_TEST(HighOrderProjection2D);
    DO_TEST(NonlinReducedIntegration);
    DO_TEST(SpectralPoisson1D);
    DO_TEST(SpectralSqrt);

    Out::root() 
      << "==================================================================="
      << endl
      << "==================================================================="
      << endl;

    if (numFail != 0)
    {
      Out::root() << "FAILURES detected!" << endl;
      Out::root() << numFail << " out of " << numFail + numPass 
                  << " tests failed" << endl;
    }
    else
    {
      Out::root() << "All tests PASSED!" << endl;
    }
    Out::root()
      << "==================================================================="
      << endl
      << "==================================================================="
      << endl;

  }
	catch(std::exception& e)
  {
    allPass = false;
    std::cerr << e.what() << std::endl;
  }
  Sundance::finalize();
  if (allPass) return 0;
  return -1;
}
