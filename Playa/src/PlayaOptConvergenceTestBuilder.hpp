#ifndef PLAYA_OPT_CONVERGENCE_TEST_BUILDER_H
#define PLAYA_OPT_CONVERGENCE_TEST_BUILDER_H


#include "PlayaOptConvergenceTestBase.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos
{
class ParameterList;
}

namespace Playa
{
using Teuchos::ParameterList;
using Teuchos::RCP;

/** 
 * Builder class for optimizer convergence tests
 *
 * @author Kevin Long
 */
class OptConvergenceTestBuilder
{
public:
  /** Construct a convergence test object from a parameter specification */
  static RCP<OptConvergenceTestBase> createConvTest(
    const ParameterList& params,
    int verb=0);
};
}

#endif
