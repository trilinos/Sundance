/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsityPattern.hpp"
#include "SundanceOut.hpp"



using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;


SumExpr::SumExpr(const RefCountPtr<ScalarExpr>& left,
                 const RefCountPtr<ScalarExpr>& right, int sign)
	: BinaryExpr(left, right, sign)
{}

bool SumExpr::isHungryDiffOp() const
{
  return leftScalar()->isHungryDiffOp() || rightScalar()->isHungryDiffOp();
}


const string& SumExpr::xmlTag() const 
{
	static string plusStr = "Plus";
	static string minusStr = "Minus";
	if (sign() < 0) return minusStr;
	return plusStr;
}

const string& SumExpr::opChar() const 
{
	static string plusStr = "+";
	static string minusStr = "-";
	if (sign() < 0) return minusStr;
	return plusStr;
}



bool SumExpr::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  TimeMonitor t(nonzeroDerivCheckTimer());
  hasNonzeroDerivCalls()++;
  /* check to see if we've already processed this node in the tree */
  if (derivHasBeenCached(d))
    {
      nonzeroDerivCacheHits()++;
      return getCachedDerivNonzeroness(d);
    }
  TimeMonitor t2(uncachedNonzeroDerivCheckTimer());
  
  /* the sum has a nonzero derivative if either operand 
   * has a nonzero derivative. */
  bool rtn = leftEvaluatable()->hasNonzeroDeriv(d) 
    || rightEvaluatable()->hasNonzeroDeriv(d);

  addDerivToCache(d, rtn);

  return rtn;
}

Array<DerivSet> SumExpr::derivsRequiredFromOperands(const DerivSet& d) const
{
  Tabs tabs;

  if (verbosity() > 1)
    {
      cerr << tabs << "SumExpr::derivsRequiredFromOperands()" << endl;
    }
  
  DerivSet leftRtn;
  DerivSet rightRtn;

  for (DerivSet::const_iterator i=d.begin(); i != d.end(); i++)
    {
      const MultipleDeriv& di = *i;
      if (leftEvaluatable()->hasNonzeroDeriv(di))
        {
          leftRtn.put(di);
        }
      if (rightEvaluatable()->hasNonzeroDeriv(di))
        {
          rightRtn.put(di);
        }
    }

  return tuple(leftRtn, rightRtn);
}

bool SumExpr::allTermsHaveTestFunctions() const
{
  return leftEvaluatable()->allTermsHaveTestFunctions()
    && rightEvaluatable()->allTermsHaveTestFunctions();
}
