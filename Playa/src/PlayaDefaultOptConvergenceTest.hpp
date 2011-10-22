#ifndef PLAYA_OPT_DEFAULT_CONVERGENCE_TEST_H
#define PLAYA_OPT_DEFAULT_CONVERGENCE_TEST_H


#include "PlayaOptConvergenceTestBase.hpp"


namespace Playa
{

/** 
 * A simple convergence test 
 *
 * @author Kevin Long
 */
class DefaultOptConvergenceTest : public OptConvergenceTestBase
{
public:
  /** */
  DefaultOptConvergenceTest(const ParameterList& params);

  /** */
  OptStatus test(const OptState& state) const ;


  /** */
  void print(std::ostream& os) const ;

private:
  int minIters_;
  int maxIters_;
  int requiredPasses_;
  double objTol_;
  double gradTol_;
  double stepTol_;
  double xTyp_;
  double fTyp_;
};

}

#endif
