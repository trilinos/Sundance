/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSumExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"



using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;


SumExpr::SumExpr(const RefCountPtr<ScalarExpr>& left,
                 const RefCountPtr<ScalarExpr>& right, int sign)
	: BinaryExpr(left, right, sign)
{
  typedef Set<int>::const_iterator setIter;

  if (isEvaluatable(left.get()) && isEvaluatable(right.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          int lod = leftEvaluatable()->orderOfSpatialDependency(d);
          int rod = rightEvaluatable()->orderOfSpatialDependency(d);
          if (lod < 0 || rod < 0) setOrderOfDependency(d, -1);
          else setOrderOfDependency(d, max(lod, rod));
        }

      Set<int> tmp = leftEvaluatable()->funcIDSet();
      tmp.merge(rightEvaluatable()->funcIDSet());
      setFuncIDSet(tmp);
      
      for (setIter i=funcIDSet().begin(); i != funcIDSet().end(); i++)
        {
          int lod = leftEvaluatable()->orderOfFunctionalDependency(*i);
          int rod = rightEvaluatable()->orderOfFunctionalDependency(*i);
          if (lod < 0 || rod < 0) setOrderOfFunctionalDependency(*i, -1);
          else setOrderOfFunctionalDependency(*i, max(lod, rod));
        }
    }
}

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


bool SumExpr::allTermsHaveTestFunctions() const
{
  return leftEvaluatable()->allTermsHaveTestFunctions()
    && rightEvaluatable()->allTermsHaveTestFunctions();
}
