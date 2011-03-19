#include "PlayaLineSearchBase.hpp"


namespace Playa
{
using Teuchos::ParameterList;

LineSearchBase::LineSearchBase(const ParameterList& params)
  : params_(params)
{
  maxSteps_ = getParameter<int>(params_, "Max Num Steps");
  minStepSize_ = getParameter<double>(params_, "Min Step Size");
  int verb = getParameter<int>(params_, "Verbosity");
  setVerb(verb);
}

}
