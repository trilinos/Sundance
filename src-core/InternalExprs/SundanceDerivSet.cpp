/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDerivSet.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

DerivSet::DerivSet()
  : SundanceUtils::Set<MultipleDeriv, increasingOrder<MultipleDeriv> >()
{}


