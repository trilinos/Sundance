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
    mediator_()
{}



void EvalManager::evalCoordExpr(const CoordExpr* expr,
                                RefCountPtr<EvalVector>& result) const 
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCoordExpr");
  mediator()->evalCoordExpr(expr, result);
}


void EvalManager::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                RefCountPtr<EvalVector>& result) const 
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalCellDiameterExpr");
  mediator()->evalCellDiameterExpr(expr, result);
}


void EvalManager::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                          const Array<MultiIndex>& mi,
                                          Array<RefCountPtr<EvalVector> >& result) const 
{
  TEST_FOR_EXCEPTION(mediator() == 0, InternalError,
                     "uninitialized mediator in "
                     "EvalManager::evalDiscreteFuncElement");

  mediator()->evalDiscreteFuncElement(expr, mi, result);
}


RefCountPtr<EvalVector> EvalManager::popVector() const
{
  return stack().popVector();
}

TempStack& EvalManager::stack()
{
  static TempStack rtn;
  return rtn;
}
