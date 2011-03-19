#include "PlayaSteepestDescent.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLineSearchBuilder.hpp"
#include "PlayaOptConvergenceTestBuilder.hpp"


namespace Playa
{
using std::endl;

SteepestDescent::SteepestDescent(
  const ParameterList& params
  )
  : LineSearchBasedOptBase(params)
{}


RCP<DirectionGeneratorBase> 
SteepestDescent::makeDirectionGenerator() const 
{
  return rcp(new SteepestDescentDirection());
}



bool SteepestDescentDirection::generateDirection(
  const RCP<ObjectiveBase>& obj,
  const Vector<double>& xCur,
  const Vector<double>& gradCur,
  const double& fCur,
  Vector<double>& p)
{
  double hInvScale = obj->getInvHScale();
  p = -hInvScale * gradCur;
  return true;
}


}
