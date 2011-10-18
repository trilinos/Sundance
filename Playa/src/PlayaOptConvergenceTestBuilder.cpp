#include "PlayaOptConvergenceTestBuilder.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"
#include "PlayaDefaultOptConvergenceTest.hpp"

namespace Playa
{

RCP<OptConvergenceTestBase> 
OptConvergenceTestBuilder::createConvTest(const ParameterList& params,
  int verb)
{
  Tabs tab(0);
  PLAYA_MSG1(verb, tab << "OptConvergenceTestBuilder::createConvTest()");
  Tabs tab1;
  PLAYA_MSG2(verb, tab1 << "params=" << params);
  
  TEUCHOS_TEST_FOR_EXCEPTION(params.name() != "Convergence Test",
    std::runtime_error, 
    "OptConvTestBuilder::createConvTest() expected parameter list named "
    "\"Convergence Test\", got name [" << params.name() << "]");

  const std::string& ctType = getParameter<string>(params, "Type");

  RCP<OptConvergenceTestBase> ct;

  if (ctType=="Default")
  {
    PLAYA_MSG2(verb, tab1 << "found Default convergence test");
    ct = rcp(new DefaultOptConvergenceTest(params));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ct.get()==0, 
    std::runtime_error, 
    "OptConvTestBuilder::createConvTest() could not construct a valid "
    "convergence test object from parameter list " << params);
    
  return ct;
}

}
