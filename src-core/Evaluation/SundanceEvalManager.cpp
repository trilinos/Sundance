/* @HEADER@ */
/* @HEADER@ */


#include "SundanceEvalManager.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceMultiIndex.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

EvalManager::EvalManager(const RefCountPtr<AbstractEvalMediator>& mediator)
  : numericalEval_(true),
    region_(),
    mediator_(mediator),
    stack_(mediator->vecSize())
{}

EvalManager::EvalManager()
  : numericalEval_(false),
    region_(),
    mediator_(),
    stack_()
{}



void EvalManager::evalCoordExpr(const CoordExpr* expr,
                                RefCountPtr<EvalVector> const & result) const 
{
  if (numericalEval_)
    {
      TEST_FOR_EXCEPTION(!result->numerical(), InternalError,
                         "EvalManager::evalCoordExpr expecting a "
                         "numerical-valued vector");
      mediator()->evalCoordExpr(expr, result.get());
    }
  else
    {
      TEST_FOR_EXCEPTION(result->numerical(), InternalError,
                         "EvalManager::evalCoordExpr expecting a "
                         "string-valued vector");
      result->setToStringValue(expr->toString());
    }
}


void EvalManager::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                       const MultiIndex& mi,
                                       RefCountPtr<EvalVector> const & result) const 
{
  if (numericalEval_)
    {
      TEST_FOR_EXCEPTION(!result->numerical(), InternalError,
                         "EvalManager::evalDiscreteFuncElement expecting a "
                         "numerical-valued vector");
      mediator()->evalDiscreteFuncElement(expr, mi, result.get());
    }
  else
    {
      TEST_FOR_EXCEPTION(result->numerical(), InternalError,
                         "EvalManager::evalDiscreteFuncElement expecting a "
                         "string-valued vector");
      if (mi.order()==0)
        {
          result->setToStringValue(expr->toString());
        }
      else
        {
          string str = "D[" + expr->toString() + ", " + mi.coordForm() + "]";
          result->setToStringValue(str);
        }
    }
}


