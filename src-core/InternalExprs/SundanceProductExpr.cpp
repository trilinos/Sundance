/* @HEADER@ */
/* @HEADER@ */


#include "SundanceProductExpr.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


ProductExpr::ProductExpr(const RefCountPtr<ScalarExpr>& left,
                         const RefCountPtr<ScalarExpr>& right)
	: BinaryExpr(left, right, 1)
{}


bool ProductExpr::isHungryDiffOp() const
{
  return rightScalar()->isHungryDiffOp();
}



const string& ProductExpr::xmlTag() const 
{
	static string timesStr = "Times";
	static string divideStr = "Divide";
	if (sign() < 0) return divideStr;
	return timesStr;
}

const string& ProductExpr::opChar() const 
{
	static string timesStr = "*";
	static string divideStr = "/";
	if (sign() < 0) return divideStr;
	return timesStr;
}



bool ProductExpr::hasNonzeroDeriv(const MultipleDeriv& d) const
{

  TimeMonitor t(nonzeroDerivCheckTimer());
  hasNonzeroDerivCalls()++;

  if (derivHasBeenCached(d))
    {
      nonzeroDerivCacheHits()++;
      return getCachedDerivNonzeroness(d);
    }

  TimeMonitor t2(uncachedNonzeroDerivCheckTimer());

  bool rtn = false;

  if (d.order()==0)
    {
      rtn = true;
    }
  else
    {
      /* If at least one term in the product rule expansion is nonzero, 
       * we have a nonzero derivative */
      Array<MultipleDeriv> leftOps;
      Array<MultipleDeriv> rightOps;

      d.productRulePermutations(leftOps, rightOps);

      for (int i=0; i<leftOps.size(); i++)
        {
          if (leftEvaluatable()->hasNonzeroDeriv(leftOps[i]) 
              && rightEvaluatable()->hasNonzeroDeriv(rightOps[i]))
            {
              rtn = true;
              break;
            }
        }
    }

  addDerivToCache(d, rtn);

  return rtn;
}

Array<DerivSet> 
ProductExpr::derivsRequiredFromOperands(const DerivSet& d) const
{
  Tabs tabs;

  if (verbosity() > 1)
    {
      cerr << tabs << "ProductExpr::derivsRequiredFromOperands()" << endl;
    }
  
  DerivSet leftRtn;
  DerivSet rightRtn;

  for (DerivSet::const_iterator i=d.begin(); i != d.end(); i++)
    {
      const MultipleDeriv& di = *i;
      if (di.order()==0)
        {
          leftRtn.put(di);
          rightRtn.put(di);
        }
      else
        {
          Array<MultipleDeriv> leftOps;
          Array<MultipleDeriv> rightOps;
          
          di.productRulePermutations(leftOps, rightOps);
          
          for (int j=0; j<leftOps.size(); j++)
            {
              if (leftEvaluatable()->hasNonzeroDeriv(leftOps[j]) 
                  && rightEvaluatable()->hasNonzeroDeriv(rightOps[j]))
                {
                  leftRtn.put(leftOps[j]);
                  rightRtn.put(rightOps[j]);
                }
            }
        }
    }

  return tuple(leftRtn, rightRtn);
}
