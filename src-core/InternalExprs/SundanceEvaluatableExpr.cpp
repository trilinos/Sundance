/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace TSFExtended;
using namespace Internal;
using namespace Internal;




EvaluatableExpr::EvaluatableExpr()
	: ScalarExpr(), 
    evaluators_(),
    sparsity_(),
    orderOfDependency_(MultiIndex::maxDim(), -1),
    nodesHaveBeenCounted_(false)
{}




RefCountPtr<SparsitySubset> 
EvaluatableExpr::sparsitySubset(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices) const 
{
  RefCountPtr<SparsitySuperset> super = sparsitySuperset(context);

  if (!super->hasSubset(multiIndices))
    {
      super->addSubset(multiIndices);
    }
  return super->subset(multiIndices);
}




RefCountPtr<SparsitySuperset> 
EvaluatableExpr::sparsitySuperset(const EvalContext& context) const 
{
  RefCountPtr<SparsitySuperset> rtn;

  if (sparsity_.containsKey(context))
    {
      rtn = sparsity_.get(context);
    }
  else
    {
      rtn = rcp(new SparsitySuperset());
      sparsity_.put(context, rtn);
    }
  return rtn;
}




const RefCountPtr<Evaluator>&
EvaluatableExpr::evaluator(const EvalContext& context) const 
{
  TEST_FOR_EXCEPTION(!evaluators_.containsKey(context), RuntimeError, 
                     "Evaluator not found for context " << context);
  return evaluators_.get(context);
}

bool EvaluatableExpr::nonzerosAreKnown(const EvalContext& context,
                                       const Set<MultiIndex>& multiIndices,
                                       bool regardFuncsAsConstant) const 
{
  NonzeroSpecifier spec(context, multiIndices, regardFuncsAsConstant);
  return knownNonzeros_.contains(spec);
}

void EvaluatableExpr::addKnownNonzero(const EvalContext& context,
                                      const Set<MultiIndex>& multiIndices,
                                      bool regardFuncsAsConstant) const 
{
  NonzeroSpecifier spec(context, multiIndices, regardFuncsAsConstant);
  knownNonzeros_.put(spec);
}



void EvaluatableExpr::evaluate(const EvalManager& mgr,
                               Array<double>& constantResults,
                               Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  evaluator(mgr.getRegion())->eval(mgr, constantResults, vectorResults);
}

void EvaluatableExpr::setupEval(const EvalContext& context) const
{
  if (!evaluators_.containsKey(context))
    {
      RefCountPtr<Evaluator> eval = rcp(createEvaluator(this, context));
      evaluators_.put(context, eval);
    }
}

void EvaluatableExpr::showSparsity(ostream& os, 
                                   const EvalContext& context) const
{
  Tabs tab0;
  os << tab0 << "Node: " << toString() << endl;
  sparsitySuperset(context)->displayAll(os);
}



int EvaluatableExpr::maxOrder(const Set<MultiIndex>& m) const 
{
  int rtn = 0;
  for (Set<MultiIndex>::const_iterator i=m.begin(); i != m.end(); i++)
    {
      rtn = max(i->order(), rtn);
    }
  return rtn;
}


const EvaluatableExpr* EvaluatableExpr::getEvalExpr(const Expr& expr)
{
  const EvaluatableExpr* rtn 
    = dynamic_cast<const EvaluatableExpr*>(expr[0].ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError,
                     "cast of " << expr 
                     << " failed in EvaluatableExpr::getEvalExpr()");
  TEST_FOR_EXCEPTION(expr.size() != 1, InternalError,
                     "non-scalar expression " << expr
                     << " in EvaluatableExpr::getEvalExpr()");

  return rtn;
}


bool EvaluatableExpr::isEvaluatable(const ExprBase* expr) 
{
  return dynamic_cast<const EvaluatableExpr*>(expr) != 0;
}


int EvaluatableExpr::countNodes() const
{
  nodesHaveBeenCounted_ = true;
  return 1;
}


