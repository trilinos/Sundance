/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_DEFS_H
#define PLAYA_DEFS_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"

#ifndef __cplusplus
#define __cplusplus
#endif

/** \brief Playa is a collection of high-level objects for linear algebra */
namespace Playa {}

/** \brief Templated overloaded operators for vectors and linear operators  */
namespace PlayaExprTemplates {}

/** \brief Functors for transformation and reduction operations */
namespace PlayaFunctors {}

namespace NOX
{
/** \brief NOXPlaya contains adapters for use of Playa vectors
* and operators in NOX solvers */
namespace NOXPlaya {}
}

#include "PlayaConfig.h"









#endif
