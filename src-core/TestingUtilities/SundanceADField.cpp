/* @HEADER@ */
/* @HEADER@ */


#include "SundanceADField.hpp"
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

ADField::ADField(const ADBasis& basis, const double& coeff)
  : basis_(basis), coeff_()
{
  double* x = new double;
  *x = coeff;
  coeff_ = rcp(x);
}

ADReal ADField::evaluate() const
{
  return *coeff_ * basis_.evaluate(evalPoint());
}

void ADField::setCoeff(const double& c)
{
  *coeff_ = c;
}

ADReal ADField::operator+(const ADReal& x) const
{
  return evaluate() + x;
}

ADReal ADField::operator+(const ADField& x) const
{
  return evaluate() + x.evaluate();
}

ADReal ADField::operator-(const ADField& x) const
{
  return evaluate() - x.evaluate();
}

ADReal ADField::operator*(const ADField& x) const
{
  return evaluate() * x.evaluate();
}


ADReal ADField::operator-(const ADReal& x) const
{
  return evaluate() - x;
}

ADReal ADField::operator+(const double& x) const
{
  return evaluate() + x;
}


ADReal ADField::operator-(const double& x) const
{
  return evaluate() - x;
}

ADReal ADField::operator-() const
{
  return -evaluate();
}

ADReal ADField::operator*(const ADReal& x) const
{
  return evaluate() * x;
}


ADReal ADField::operator*(const double& x) const
{
  return evaluate() * x;
}

