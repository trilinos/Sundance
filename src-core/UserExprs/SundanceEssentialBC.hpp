/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ESSENTIALBC_H
#define SUNDANCE_ESSENTIALBC_H

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
  Expr EssentialBC(const Handle<CellFilterBase>& domain,
                   const Expr& integrand);

  /** \relates Expr */
  Expr EssentialBC(const Handle<CellFilterBase>& domain,
                   const Expr& integrand,
                   const Handle<QuadratureFamilyBase>& quad);

  
}

#endif
