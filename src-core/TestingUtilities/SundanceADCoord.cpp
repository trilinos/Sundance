/* @HEADER@ */
/* @HEADER@ */


#include "SundanceADCoord.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceTesting;
using namespace TSFExtended;
using namespace Teuchos;
using namespace std;

ADCoord::ADCoord(int dir)
  : dir_(dir)
{}

ADReal ADCoord::evaluate() const
{
  Point grad(0.0, 0.0, 0.0);
  grad[dir_]=1.0;
  ADReal rtn = ADReal(ADField::evalPoint()[dir_], grad);
  return rtn;
}



ADReal ADCoord::operator+(const ADReal& x) const
{
  return evaluate() + x;
}

ADReal ADCoord::operator+(const ADCoord& x) const
{
  return evaluate() + x.evaluate();
}

ADReal ADCoord::operator+(const ADField& x) const
{
  return evaluate() + x.evaluate();
}

ADReal ADCoord::operator+(const double& x) const
{
  return evaluate() + x;
}



ADReal ADCoord::operator-(const ADCoord& x) const
{
  return evaluate() - x.evaluate();
}

ADReal ADCoord::operator-(const ADReal& x) const
{
  return evaluate() - x;
}

ADReal ADCoord::operator-(const ADField& x) const
{
  return evaluate() - x;
}

ADReal ADCoord::operator-(const double& x) const
{
  return evaluate() - x;
}



ADReal ADCoord::operator-() const
{
  return -evaluate();
}

ADReal ADCoord::operator*(const ADCoord& x) const
{
  return evaluate() * x.evaluate();
}

ADReal ADCoord::operator*(const ADReal& x) const
{
  return evaluate() * x;
}

ADReal ADCoord::operator*(const double& x) const
{
  return evaluate() * x;
}

ADReal ADCoord::operator*(const ADField& x) const
{
  return evaluate() * x.evaluate();
}

