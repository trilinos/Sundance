/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDeriv.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
 

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


bool Deriv::operator<(const Deriv& other) const 
{
  return ptr()->lessThan(other);
}

bool Deriv::operator==(const Deriv& other) const
{
  if (*this < other || other < *this) return false;
}


string Deriv::toString() const 
{
  return ptr()->toString();
}

const FunctionalDeriv* Deriv::funcDeriv() const 
{
  return dynamic_cast<FunctionalDeriv*>(ptr().get());
}

const CoordDeriv* Deriv::coordDeriv() const 
{
  return dynamic_cast<CoordDeriv*>(ptr().get());
}

bool Deriv::isTestFunction() const 
{
  if (!isFunctionalDeriv()) return false;

  const FunctionalDeriv* f = funcDeriv();

  TEST_FOR_EXCEPTION(f==0, RuntimeError,
                     "contents not a functional derivative in Deriv::isTestFunction()");

  const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(f->func());
  return t != 0;
}

bool Deriv::isUnknownFunction() const 
{
  if (!isFunctionalDeriv()) return false;

  const FunctionalDeriv* f = funcDeriv();

  TEST_FOR_EXCEPTION(f==0, RuntimeError,
                     "contents not a functional derivative in Deriv::isTestFunction()");

  const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(f->func());
  return u != 0;
}
