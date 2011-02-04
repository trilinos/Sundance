/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaDebug.hpp"

namespace Playa
{
#ifdef TEUCHOS_DEBUG
bool Debug::on = true;
#else
bool Debug::on = false;
#endif
}

