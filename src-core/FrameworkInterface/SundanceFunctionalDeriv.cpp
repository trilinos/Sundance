/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceOrderedTuple.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

FunctionalDeriv::FunctionalDeriv(const FuncElementBase* func,
                                 const MultiIndex& mi)
  : DerivBase(), func_(func), mi_(mi)
{;}

bool FunctionalDeriv::lessThan(const Deriv& other) const
{
  /* First compare types: spatial derivs are ranked lower than 
   * functional derivs, so if the other guy is a spatial deriv, I'm higher */
  if (other.isCoordDeriv()) return false;

  /* We can now safely cast to a functional deriv and compare contents */

  const FunctionalDeriv* f = other.funcDeriv();

  if (funcID() < f->funcID()) return true;
  if (funcID() > f->funcID()) return false;
  if (multiIndex() < f->multiIndex()) return true;
  return false;
}

Deriv FunctionalDeriv::derivWrtMultiIndex(const MultiIndex& mi) const
{
  return new FunctionalDeriv(func_, mi_ + mi);
}

string FunctionalDeriv::toString() const 
{
  if (mi_.order()==0)
    {
      return func_->name();
    }
  return "D[" + func_->name() + ", " + mi_.coordForm() + "]";
}


