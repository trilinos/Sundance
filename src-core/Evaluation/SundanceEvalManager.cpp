/* @HEADER@ */
/* @HEADER@ */


#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
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

EvalManager::EvalManager()
  : region_(),
    mediator_(),
    stack_()
{}



void EvalManager::evalCoordExpr(const CoordExpr* expr,
                                RefCountPtr<EvalVector> const & result) const 
{
  if (mediator() != 0)
    {
      result->setToVectorValue();
      mediator()->evalCoordExpr(expr, result.get());
    }
  else
    {
      result->setToVectorValue();
      result->resize(1);
      result->setElement(0, 1.0);
    }
  
  if (result->shadowOps())
    {
      result->setStringValue(expr->toString());
    }
}


void EvalManager::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                       const MultiIndex& mi,
                                       RefCountPtr<EvalVector> const & result) const 
{
  if (mediator() != 0)
    {
      result->setToVectorValue();
      mediator()->evalDiscreteFuncElement(expr, mi, result.get());
    }
  else
    {
      result->setToVectorValue();
      result->resize(1);
      result->setElement(0, 1.0);
    }
  
  if (result->shadowOps())
    {
      if (mi.order()==0)
        {
          result->setStringValue(expr->toString());
        }
      else
        {
          string str = "D[" + expr->toString() + ", " + mi.coordForm() + "]";
          result->setStringValue(str);
        }
    }
}


