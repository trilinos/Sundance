/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ESSENTIALBC_H
#define SUNDANCE_ESSENTIALBC_H

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
  Expr EssentialBC(const Handle<CellFilterStub>& domain,
                   const Expr& integrand);

  /** \relates Expr */
  Expr EssentialBC(const Handle<CellFilterStub>& domain,
                   const Expr& integrand,
                   const Handle<QuadratureFamilyStub>& quad);

  
}

#endif
