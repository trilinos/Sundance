#include "PlayaLineSearchBuilder.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"
#include "PlayaSimpleBacktracking.hpp"

namespace Playa
{

RCP<LineSearchBase> 
LineSearchBuilder::createLineSearch(const ParameterList& params,
  int verb)
{
  Tabs tab(0);
  PLAYA_MSG1(verb, tab << "LineSearchBuilder::createLineSearch()");
  Tabs tab1;
  PLAYA_MSG2(verb, tab1 << "params=" << params);
  
  TEUCHOS_TEST_FOR_EXCEPTION(params.name() != "Line Search",
    std::runtime_error, 
    "LineSearchBuilder::getLineSearch() expected parameter list named "
    "\"Line Search\", got name [" << params.name() << "]");

  const std::string& lsType = getParameter<string>(params, "Type");

  RCP<LineSearchBase> ls;

  if (lsType=="Simple Backtracking")
  {
    PLAYA_MSG2(verb, tab1 << "found Simple Backtracking LS");
    ls = rcp(new SimpleBacktracking(params));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ls.get()==0, 
    std::runtime_error, 
    "LineSearchBuilder::getLineSearch() could not construct a valid line "
    "search object from parameter list " << params);
    
  return ls;
}

}
