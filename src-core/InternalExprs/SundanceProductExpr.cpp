/* @HEADER@ */
/* @HEADER@ */


#include "SundanceProductExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


ProductExpr::ProductExpr(const RefCountPtr<ScalarExpr>& left,
                         const RefCountPtr<ScalarExpr>& right)
	: BinaryExpr(left, right, 1)
{
  if (isEvaluatable(left.get()) && isEvaluatable(right.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          int lod = leftEvaluatable()->orderOfDependency(d);
          int rod = rightEvaluatable()->orderOfDependency(d);
          if (lod < 0 || rod < 0) setOrderOfDependency(d, -1);
          else setOrderOfDependency(d, lod+rod);
        }
    }
}


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

void ProductExpr::findChildMultiIndexSets(const Set<MultiIndex>& miSet,
                                          Set<MultiIndex>& miLeft,
                                          Set<MultiIndex>& miRight) const
{
  for (Set<MultiIndex>::const_iterator i=miSet.begin(); 
       i != miSet.end(); i++)
    {
      const MultiIndex& mi = *i;

      TEST_FOR_EXCEPTION(mi.order() > 1, RuntimeError,
                         "multiindex order restricted to 0,1");

      /* if first derivatives are needed, by the product rule we
       * may need zero-order derivs also */
      if (mi.order() == 1) 
        {
          int d = mi.firstOrderDirection();
          if (leftEvaluatable()->orderOfDependency(d) != 0)
            {
              miRight.put(MultiIndex());
              miLeft.put(mi);
            }
          if (rightEvaluatable()->orderOfDependency(d) != 0)
            {
              miLeft.put(MultiIndex());
              miRight.put(mi);
            }
        }
      else
        {
          miLeft.put(mi);
          miRight.put(mi);
        }
    }
}





void ProductExpr::findNonzeros(const EvalContext& context,
                               const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                               bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for product expr" 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  Set<MultiIndex> miLeft;
  Set<MultiIndex> miRight;
  findChildMultiIndexSets(multiIndices, miLeft, miRight);

  Set<MultiSet<int> > childFuncIDs = findChildFuncIDSet(activeFuncIDs,
                                                        allFuncIDs);

  int maxSpatialOrder = maxOrder(multiIndices);
  int maxDiffOrder = context.topLevelDiffOrder() + maxSpatialOrder;



  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: getting left operand's nonzeros");
  leftEvaluatable()->findNonzeros(context, miLeft,
                                  childFuncIDs,
                                  allFuncIDs,
                                  regardFuncsAsConstant);

  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: getting right operand's nonzeros");
  rightEvaluatable()->findNonzeros(context, miRight,
                                   childFuncIDs,
                                   allFuncIDs,
                                   regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> leftSparsity 
    = leftEvaluatable()->sparsitySubset(context, miLeft);

  RefCountPtr<SparsitySubset> rightSparsity 
    = rightEvaluatable()->sparsitySubset(context, miRight);

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  for (int i=0; i<leftSparsity->numDerivs(); i++)
    {
      const MultipleDeriv& dLeft = leftSparsity->deriv(i);
     
      for (int j=0; j<rightSparsity->numDerivs(); j++)
        {
          const MultipleDeriv& dRight = rightSparsity->deriv(j);

          /* Skip combinations of functional derivatives that contribute
           * only to derivatives of an order we don't need */
          if (dRight.order() + dLeft.order() > maxDiffOrder) continue;

          /* Skip combinations of spatial derivs of greater order
           * than the max order of our multiindices */
          if (dRight.spatialOrder() + dLeft.spatialOrder() > maxSpatialOrder) continue;

          /*
           * Skip combinations that do not contribute to the
           * variational derivatives required at this point.
           */
          MultiSet<int> funcs;
          for (MultipleDeriv::const_iterator k=dLeft.begin();
               k != dLeft.end(); k++)
            {
              const Deriv& d = *k;
              if (d.isFunctionalDeriv()) 
                {
                  int fid = d.funcDeriv()->funcID();
                  funcs.put(fid);
                }
            }
          for (MultipleDeriv::const_iterator k=dRight.begin();
               k != dRight.end(); k++)
            {
              const Deriv& d = *k;
              if (d.isFunctionalDeriv()) 
                {
                  int fid = d.funcDeriv()->funcID();
                  funcs.put(fid);
                }
            }
          
          if (!activeFuncIDs.contains(funcs)) continue;


          /* The current left and right nonzero functional derivatives
           * dLeft, dRight will contribute to the dLeft*dRight 
           * functional derivative */
          MultipleDeriv productDeriv = dLeft;
          for (MultipleDeriv::const_iterator k=dRight.begin();
               k != dRight.end(); k++)
            {
              productDeriv.insert(*k);
            }
          /* Use the more general of the two operands' states */
          DerivState newState = max(leftSparsity->state(i), 
                                    rightSparsity->state(j));
          subset->addDeriv(productDeriv, newState);
        }
    }

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant);
}
