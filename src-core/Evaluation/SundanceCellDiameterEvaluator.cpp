/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellDiameterExprEvaluator::CellDiameterExprEvaluator(const CellDiameterExpr* expr, 
                                       const EvalContext& context)
  : SubtypeEvaluator<CellDiameterExpr>(expr, context), 
    stringRep_(expr->toString())
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing cell diameter expr evaluator for " 
                    << expr->toString());
  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  TEST_FOR_EXCEPTION(sparsity()->numDerivs() > 1, InternalError,
                     "CellDiameterExprEvaluator ctor found a sparsity table "
                     "with more than one entry. The bad sparsity table is "
                     << *sparsity());

  /* 
   * There is only one possible entry in the nozeros table for a
   * cell diameter expression: a zeroth derivative.
   */
  
  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity()->deriv(i);

      TEST_FOR_EXCEPTION(d.order()!=0, InternalError,
                         "CellDiameterExprEvaluator ctor found an entry in the "
                         "sparsity superset that is not a zeroth-order derivative. "
                         "The bad entry is " << sparsity()->deriv(i) 
                         << ". The superset is " 
                         << *sparsity());
      addVectorIndex(i, 0);
    }
}



void CellDiameterExprEvaluator::internalEval(const EvalManager& mgr,
                                             Array<double>& constantResults,
                                             Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(cellDiameterEvalTimer());
  Tabs tabs;

  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---CellDiameterExprEvaluator---");

  if (verbosity() > 1)
    {
      cerr << tabs << "CellDiameterExprEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      cerr << tabs << "sparsity = " << endl << *sparsity() << endl;
    }

  vectorResults.resize(1);
  vectorResults[0] = mgr.popVector();
  mgr.evalCellDiameterExpr(expr(), vectorResults[0]);
  if (EvalVector::shadowOps()) vectorResults[0]->setString(stringRep_);

  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "results " << endl;
      sparsity()->print(cerr, vectorResults,
                            constantResults);
    }
}

