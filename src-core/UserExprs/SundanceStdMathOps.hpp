/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDMATHOPS_H
#define SUNDANCE_STDMATHOPS_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  /** \relates Expr */
  Expr sin(const Expr& expr);
  
  /** \relates Expr */
  Expr cos(const Expr& expr);

  /** \relates Expr */
  Expr exp(const Expr& expr);

}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
