/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRAL_H
#define SUNDANCE_INTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceCellFilterBase.hpp"
#include "TSFHandle.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace TSFExtended;
  using namespace FrameworkInterface;
  using namespace Teuchos;

  /** \relates Expr */
  Expr Integral(const Handle<CellFilterBase>& domain,
                const Expr& integrand);

  /** \relates Expr */
  Expr Integral(const Handle<CellFilterBase>& domain,
                const Expr& integrand,
                const Handle<QuadratureFamilyBase>& quad);

  
}

#endif
