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


StdFwkEvalMediator::StdFwkEvalMediator(const Mesh& mesh, int cellDim)
  : mesh_(mesh),
    cellDim_(cellDim),
    cellType_(NullCell),
    cellLID_(),
    J_(rcp(new CellJacobianBatch())),
    cacheIsValid_(false),
    jCacheIsValid_(false)
{;}

void StdFwkEvalMediator::getJacobians(RefCountPtr<CellJacobianBatch>& J) const
{
  if (!jCacheIsValid_)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "computing new Jacobians");
      mesh().getJacobians(cellDim(), *cellLID(), *J_);
      jCacheIsValid_ = true;
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "reusing cached Jacobians");
    }
  J = J_;
}

