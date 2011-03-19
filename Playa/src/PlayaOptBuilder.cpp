#include "PlayaOptBuilder.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"
#include "PlayaSteepestDescent.hpp"
#include "PlayaBasicLMBFGS.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

namespace Playa
{

RCP<UnconstrainedOptimizerBase> 
OptBuilder::createOptimizer(const string& filename,
  int verb)
{
  ParameterXMLFileReader reader(filename);
  ParameterList params = reader.getParameters();
  return createOptimizer(params);
}

RCP<UnconstrainedOptimizerBase> 
OptBuilder::createOptimizer(const ParameterList& params,
  int verb)
{
  Tabs tab(0);
  PLAYA_MSG1(verb, tab << "OptBuilder::createOptimizer()");
  Tabs tab1;
  PLAYA_MSG2(verb, tab1 << "params=" << params);
  
  TEST_FOR_EXCEPTION(params.name() != "Optimizer",
    std::runtime_error, 
    "Optimizer::getOptimizer() expected parameter list named "
    "\"Optimizer\", got name [" << params.name() << "]");

  const std::string& optType = getParameter<string>(params, "Type");

  RCP<UnconstrainedOptimizerBase> opt;

  if (optType=="Steepest Descent")
  {
    PLAYA_MSG2(verb, tab1 << "found SteepestDescent");
    opt = rcp(new SteepestDescent(params));
  }
  else if (optType=="Basic LMBFGS")
  {
    PLAYA_MSG2(verb, tab1 << "found Basic LMBFGS");
    opt = rcp(new BasicLMBFGS(params));
  }

  TEST_FOR_EXCEPTION(opt.get()==0, 
    std::runtime_error, 
    "OptBuilder::createOptimizer() could not construct a valid "
    "optimizer object from parameter list " << params);
    
  return opt;
}

}
