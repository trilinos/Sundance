/* @HEADER@ */
/* @HEADER@ */

#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;




void StdFwkEvalMediator::setCellBatch(const RefCountPtr<Array<int> >& cellLID,
                                      RefCountPtr<CellJacobianBatch>& J) 
{
  cellLID_ = cellLID; 
  cacheIsValid() = false; 
  jCacheIsValid_=false;
  mesh_.getJacobians(cellDim(), *cellLID, *J);
  J_ = J;

  /* mark the function caches as invalid */
  Map<const DiscreteFunction*, bool>::iterator iter;
  for (iter = fCacheIsValid_.begin(); iter != fCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = dfCacheIsValid_.begin(); iter != dfCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
  for (iter = localValueCacheIsValid_.begin(); iter != localValueCacheIsValid_.end(); iter++)
    {
      iter->second = false;
    }
}
