/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCoordDeriv.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDeriv.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

CoordDeriv::CoordDeriv(int dir)
  : DerivBase(), dir_(dir)
{}

bool CoordDeriv::lessThan(const Deriv& other) const
{
  /* First compare types: spatial derivs are ranked lower than 
   * functional derivs, so if the other guy is a func deriv, I'm lower */
  if (other.isFunctionalDeriv()) return true;

  /* we can now safely cast the other guy to a spatial deriv and compare
   * directions */
  const CoordDeriv* s = other.coordDeriv();

  return dir_ < s->dir_;
}

string CoordDeriv::toString() const 
{
  return CoordExpr::coordName(dir_, "");
}


