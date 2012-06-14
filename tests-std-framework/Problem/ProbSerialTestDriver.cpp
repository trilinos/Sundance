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
      failures.append(#testName); \
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
    Array<string> failures;

    DO_TEST(NonlinearPartialDomain);
    DO_TEST(NonlinearPeriodic1D);
    DO_TEST(LinearPeriodic1D);
    DO_TEST(PoissonOnDisk);
    DO_TEST(TetQuadTransformationTest);
//    DO_TEST(DiscFunc3D);
    DO_TEST(Kepler);
    DO_TEST(CNBugTest);
    DO_TEST(EdgeDFTest);
    DO_TEST(AToCDensitySample);
    DO_TEST(DuffingFloquet);
    DO_TEST(SecondOrderFloquet);

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
      Out::root() << "failures are: " << failures << endl;
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
