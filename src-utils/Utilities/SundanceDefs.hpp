/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DEFS_H
#define SUNDANCE_DEFS_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"

#ifndef __cplusplus
#define __cplusplus
#endif

#ifdef HAVE_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package 
 * and need to be undef'd here to avoid warnings when this file is 
 * included from another package.
 * KL 11/25/02
 */   
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "SundanceConfig.h"

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif


#endif /* HAVE_CONFIG_H */

namespace Teuchos {;}

namespace SundanceUtils {;}

namespace SundanceCore 
{
  namespace Internal{;}
}

namespace SundanceStdMesh
{
  namespace Internal{;}
}

#endif
