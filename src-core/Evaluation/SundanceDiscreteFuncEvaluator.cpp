/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


DiscreteFuncElementEvaluator
::DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<DiscreteFuncElement>(expr, context), 
    mi_(sparsity()->numDerivs()),
    miToIndexMap_(),
    stringReps_()
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "initializing discrete func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  string funcName = expr->name();

  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      /* Make sure that every derivative we're to evaluate is either
      * zero-order or a spatial derivative */
      if (sparsity()->deriv(i).order()==0) 
        {
          mi_[i] = MultiIndex();
        }
      else 
        {
      
          TEST_FOR_EXCEPTION(!sparsity()->isSpatialDeriv(i), InternalError,
                             "DiscreteFuncElementEvaluator ctor found "
                             "an entry in the sparsity superset that is not "
                             "a spatial derivative. "
                             "The bad entry is " << sparsity()->deriv(i) 
                             << ". The superset is " 
                             << *sparsity());

          mi_[i] = sparsity()->multiIndex(i);
        }
      addVectorIndex(i,i);
      TEST_FOR_EXCEPTION(miToIndexMap_.containsKey(mi_[i]), InternalError,
                         "DiscreteFuncElementEvaluator ctor detected a "
                         "duplicate multiindex");

      miToIndexMap_.put(mi_[i], i);

      if (mi_[i].order()==0)
        {
          stringReps_.append(funcName);
        }
      else
        {
          int dir = mi_[i].firstOrderDirection();
          string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
    }

  
}

bool DiscreteFuncElementEvaluator::hasMultiIndex(const MultiIndex& mi) const
{
  return miToIndexMap_.containsKey(mi);
}

int DiscreteFuncElementEvaluator::miIndex(const MultiIndex& mi) const
{
  return miToIndexMap_.get(mi);
}


void DiscreteFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(discreteFuncEvalTimer());
  Tabs tabs;

  if (verbosity() > VerbSilent)
    {
      cerr << tabs << "DiscreteFuncElementEvaluator::eval: expr=" << expr()->toString() 
           << endl;
    }

  vectorResults.resize(mi_.size());
  for (unsigned int i=0; i<mi_.size(); i++)
    {
      vectorResults[i] = mgr.popVector();
      TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                         InternalError,
                         "invalid evaluation vector allocated in "
                         "DiscreteFuncElementEvaluator::internalEval()");
      vectorResults[i]->setString(stringReps_[i]);
    }
  mgr.evalDiscreteFuncElement(expr(), mi_, vectorResults);
  
  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "results " << endl;
      sparsity()->print(cerr, vectorResults,
                            constantResults);
    }
}

