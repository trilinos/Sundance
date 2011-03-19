#ifndef PLAYA_OPT_CONVERGENCE_TEST_BASE_H
#define PLAYA_OPT_CONVERGENCE_TEST_BASE_H


#include "PlayaObjectWithVerbosity.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaOptState.hpp"


namespace Playa
{
class OptState;
using Teuchos::ParameterList;

/** 
 * Base class for convergence tests for optimization algorithms
 *
 * @author Kevin Long
 */
class OptConvergenceTestBase : public ObjectWithVerbosity,
                               public Printable
{
public:
  /** */
  OptConvergenceTestBase(const ParameterList& params)
    {
      int verb = params.get<int>("Verbosity");
      setVerb(verb);
    }

  /** */
  virtual ~OptConvergenceTestBase() {}

  /** */
  virtual OptStatus test(const OptState& state) const = 0 ;
};

}

#endif
