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
  // Set<MultiSet<int> > rtn;

//   for (Set<MultiSet<int> >::const_iterator 
//          i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
//     {
//       const MultiSet<int>& ms = *i;
//       /* if any derivs of order > 0 are requested, we'll need the zero-order
//        * deriv also for any nonlinear operation. */
//       if (ms.size() > 0) rtn.put(MultiSet<int>());

      
//     }
}

