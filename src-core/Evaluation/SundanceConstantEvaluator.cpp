/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;



ConstantEvaluator::ConstantEvaluator(const SpatiallyConstantExpr* expr, 
                                     const EvalContext& context)
  : SubtypeEvaluator<SpatiallyConstantExpr>(expr, context)
{
  /*
   * There is only one possible nonzero derivative of this expression: the
   * zeroth-order derivative. 
   *
   * There's nothing to do in this ctor other than running some sanity checks.
   */

  TEST_FOR_EXCEPTION(sparsity()->numDerivs() != 1, InternalError,
                     "ConstantEvaluator ctor found a sparsity table "
                     "without exactly one entry. The bad sparsity table is "
                     << *sparsity());

  const MultipleDeriv& d = sparsity()->deriv(0);

  TEST_FOR_EXCEPTION(d.order() != 0, InternalError,
                     "ConstantEvaluator ctor found a nonzero derivative "
                     "of order greater than zero. The bad sparsity table is "
                     << *sparsity());

  addConstantIndex(0,0);
}





void ConstantEvaluator::internalEval(const EvalManager& mgr,
                                     Array<double>& constantResults,
                                     Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(constantEvalTimer());
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---ConstantEvaluator---");

  constantResults.resize(1);
  constantResults[0] = expr()->value();
}
