/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_INTEGRAL_H
#define SUNDANCE_INTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "TSFHandle.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace TSFExtended;
  using namespace Internal;
  using namespace Teuchos;

  /** \relates Expr */
  Expr Integral(const Handle<CellFilterStub>& domain,
                const Expr& integrand);

  /** \relates Expr */
  Expr Integral(const Handle<CellFilterStub>& domain,
                const Expr& integrand,
                const Handle<QuadratureFamilyStub>& quad);

  
}

#endif
