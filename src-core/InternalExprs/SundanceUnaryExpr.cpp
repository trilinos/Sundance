/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


UnaryExpr::UnaryExpr(const RefCountPtr<ScalarExpr>& arg)
	: ExprWithChildren(tuple(arg)), allActiveFuncs_()
{}


void UnaryExpr::addActiveFuncs(const EvalContext& context,
                               const Set<MultiSet<int> >& activeFuncIDs) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "Expr " << toString() << " adding active funcs "
                     << activeFuncIDs);
  if (allActiveFuncs_.containsKey(context))
    {
      allActiveFuncs_[context].merge(activeFuncIDs);        
    }
  else
    {
      allActiveFuncs_.put(context, activeFuncIDs);
    }

}

const Set<MultiSet<int> >& UnaryExpr::getActiveFuncs(const EvalContext& context) const 
{
  TEST_FOR_EXCEPTION(!allActiveFuncs_.containsKey(context), 
                     InternalError,
                     "context " << context << " does not exist in UnaryExpr::getActiveFuncs()");
  return allActiveFuncs_.get(context);
}
