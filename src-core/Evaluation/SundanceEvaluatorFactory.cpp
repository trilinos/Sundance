/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceInstructionCachingEvaluator.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace SundanceCore::Internal;


EvaluatorFactory::EvaluatorFactory()
{;}

Evaluator* EvaluatorFactory::commonCreate(const EvaluatableExpr* expr,
                                          const EvalContext& context,
                                          int topLevelDiffOrder) const
{
  const CoordExpr* c = dynamic_cast<const CoordExpr*>(expr);
  if (c != 0)
    {
      return new CoordExprEvaluator(c, context, topLevelDiffOrder);
    }
  
  const SpatiallyConstantExpr* sc 
    = dynamic_cast<const SpatiallyConstantExpr*>(expr);
  if (sc != 0)
    {
      return new ConstantEvaluator(sc, context, topLevelDiffOrder);
    }

  const SymbolicFuncElement* u
    = dynamic_cast<const SymbolicFuncElement*>(expr);
  if (u != 0)
    {
      return new SymbolicFuncElementEvaluator(u, context, topLevelDiffOrder);
    }
  
  const DiscreteFuncElement* df
    = dynamic_cast<const DiscreteFuncElement*>(expr);
  if (df != 0)
    {
      return new DiscreteFuncElementEvaluator(df, context, topLevelDiffOrder);
    }

  TEST_FOR_EXCEPTION(true, InternalError,
                     "EvaluatorFactory::commonCreate() could not create an "
                     "evaluator for " << expr->toString());

  return 0;
}


RefCountPtr<EvaluatorFactory>&  EvaluatorFactory::defaultEvaluator()
{
  static RefCountPtr<EvaluatorFactory> rtn 
    = rcp(new InstructionCachingEvaluatorFactory());
  return rtn;
}
