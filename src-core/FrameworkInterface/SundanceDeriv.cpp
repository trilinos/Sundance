/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDeriv.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceExceptions.hpp"

 

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


