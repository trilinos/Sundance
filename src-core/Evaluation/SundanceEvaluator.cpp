/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluator.hpp"
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



Evaluator::Evaluator()
  : numClients_(0),
    numCalls_(0),
    vectorResultCache_(),
    constantResultCache_(),
    constantIndexMap_(),
    vectorIndexMap_(),
    vectorIndices_(),
    constantIndices_()
{}



void Evaluator::eval(const EvalManager& mgr,
                     Array<double>& constantResults,
                     Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  Tabs tabs;

  if (numCalls_ == 0)
    {
      internalEval(mgr, constantResultCache_, vectorResultCache_);
    }
  
  numCalls_++;

  /* Go ahead and copy the constant results every time, 
   * since this is cheap */
  constantResults = constantResultCache_;

  /* If all clients have called, we can return the original data
   * which can then be changed by the client. */
  if (numCalls_ == numClients_)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "surrendering cached results");
      vectorResults = vectorResultCache_;
    }
  else /* Otherwise, make a copy of the data */
    {
      SUNDANCE_VERB_MEDIUM(tabs << "cloning cached results");
      vectorResults.resize(vectorResultCache_.size());
      for (unsigned int i=0; i < vectorResults.size(); i++)
        {
          vectorResults[i] = vectorResultCache_[i]->clone();
        }
    }
}


void Evaluator::addConstantIndex(int index, int constantIndex)
{
  TEST_FOR_EXCEPTION(constantIndexMap_.containsKey(index), InternalError,
                     "duplicate index " << index 
                     << " found in Evaluator::addConstantIndex");
  constantIndexMap_.put(index, constantIndex);
  constantIndices_.append(index);
}

void Evaluator::addVectorIndex(int index, int vectorIndex)
{
  TEST_FOR_EXCEPTION(vectorIndexMap_.containsKey(index), InternalError,
                     "duplicate index " << index 
                     << " found in Evaluator::addVectorIndex");
  vectorIndexMap_.put(index, vectorIndex);
  vectorIndices_.append(index);
}

