/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

SpatiallyConstantExpr::SpatiallyConstantExpr(const double& value)
	: EvaluatableExpr(), value_(value)
{
  for (int d=0; d<MultiIndex::maxDim(); d++) 
    {
      setOrderOfDependency(d, 0);
    }
}




void SpatiallyConstantExpr::findNonzeros(const EvalContext& context,
                                         const Set<MultiIndex>& multiIndices,
                                         bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for constant" 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  if (nonzerosAreKnown(context, multiIndices, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);
  subset->addDeriv(MultipleDeriv(), ConstantDeriv);

  addKnownNonzero(context, multiIndices, regardFuncsAsConstant);
}
