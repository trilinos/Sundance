/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNonlinearExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;

Set<MultiSet<int> > 
NonlinearExpr::findChildFuncIDSet(const Set<MultiSet<int> >& activeFuncIDs,
                                  const Set<int>& allFuncIDs) const 
{
  return activeFuncIDs;
}

