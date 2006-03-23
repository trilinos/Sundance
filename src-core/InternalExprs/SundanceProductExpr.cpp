/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */


#include "SundanceProductExpr.hpp"
#include "SundanceProductEvaluator.hpp"
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
          int lod = leftEvaluatable()->orderOfSpatialDependency(d);
          int rod = rightEvaluatable()->orderOfSpatialDependency(d);
          if (lod < 0 || rod < 0) setOrderOfDependency(d, -1);
          else setOrderOfDependency(d, lod+rod);
        }
      
      const Set<MultiSet<int> >& leftFuncs = leftEvaluatable()->funcIDSet();
      const Set<MultiSet<int> >& rightFuncs = rightEvaluatable()->funcIDSet();
      typedef Set<MultiSet<int> >::const_iterator iter;
      
      for (iter i=leftFuncs.begin(); i != leftFuncs.end(); i++)
        {
          const MultiSet<int>& fl = *i;
          for (iter j=rightFuncs.begin(); j != rightFuncs.end(); j++)
            {
              const MultiSet<int>& fr = *j;
              if (fl.size() + fr.size() > maxFuncDiffOrder()) continue;
              addFuncIDCombo(fr.merge(fl));
            }
        }
    }
}


Evaluator* ProductExpr::createEvaluator(const EvaluatableExpr* expr,
					const EvalContext& context) const
{
  return new ProductEvaluator(dynamic_cast<const ProductExpr*>(expr), context);
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
          if (leftEvaluatable()->orderOfSpatialDependency(d) != 0)
            {
              miRight.put(MultiIndex());
              miLeft.put(mi);
            }
          if (rightEvaluatable()->orderOfSpatialDependency(d) != 0)
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


void ProductExpr
::findChildActiveFuncs(const Set<MultiSet<int> >& funcIDs,
                       Set<MultiSet<int> >& leftFuncs,
                       Set<MultiSet<int> >& rightFuncs) const 
{
  Tabs tabs;

  typedef Set<MultiSet<int> >::const_iterator iter;
  
  //  SUNDANCE_VERB_HIGH(tabs << "doing product rule");

  for (iter f=funcIDs.begin(); f != funcIDs.end(); f++)
    {
      Tabs tab1;
      //  SUNDANCE_VERB_HIGH(tab1 << "partitioning " << *f);
      const MultiSet<int>& d = *f;
      int n = d.size();
      if (n==0)
        {
          Tabs tab2;
          leftFuncs.put(MultiSet<int>());
          rightFuncs.put(MultiSet<int>());
          //  SUNDANCE_VERB_HIGH(tab2 << "found L={}, R={}");
        }
      else
        {
          Tabs tab2;
          int p2 = MultipleDeriv::pow2(n);
          for (int i=0; i<p2; i++)
            {
              MultiSet<int> left;
              MultiSet<int> right;
              Array<int> bits = MultipleDeriv::bitsOfAnInteger(i, n);
              int k=0;
              for (MultiSet<int>::const_iterator 
                     j=d.begin(); j != d.end(); j++, k++)
                {
                  if (bits[k]==true)
                    {
                      left.put(*j);
                    }
                  else
                    {
                      right.put(*j);
                    }
                }
              //  SUNDANCE_VERB_HIGH(tab2 << "found L=" << left.toString()
              //                   << " R=" << right.toString());
              if (leftEvaluatable()->funcIDSet().contains(left))
                {
                  if (rightEvaluatable()->funcIDSet().contains(right)) 
                    rightFuncs.put(right);
                }
              else
                {
                  //  SUNDANCE_VERB_HIGH(tab2 << "skipping R=" << right.toString()
                  //                 << " since L=" << left.toString()
                  //                 << " is zero");
                }
              if (rightEvaluatable()->funcIDSet().contains(right))
                {
                  if (leftEvaluatable()->funcIDSet().contains(left))
                    leftFuncs.put(left);
                }
              else
                {
                  // SUNDANCE_VERB_HIGH(tab2 << "skipping L=" << left.toString()
                  //                 << " since R=" << right.toString()
                  //                 << " is zero");
                }
            }
        }
    }
}
                       
                                       





void ProductExpr::findNonzeros(const EvalContext& context,
                               const Set<MultiIndex>& multiIndices,
                               const Set<MultiSet<int> >& inputActiveFuncIDs,
                               bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for product expr" 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  Set<MultiSet<int> > activeFuncIDs = filterActiveFuncs(inputActiveFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs, false);

  Set<MultiIndex> miLeft;
  Set<MultiIndex> miRight;
  findChildMultiIndexSets(multiIndices, miLeft, miRight);

  Set<MultiSet<int> > leftFuncIDs;
  Set<MultiSet<int> > rightFuncIDs;

  findChildActiveFuncs(activeFuncIDs,
                       leftFuncIDs,
                       rightFuncIDs);

  if (miLeft.size() > 0) rightFuncIDs.put(MultiSet<int>());
  if (miRight.size() > 0) leftFuncIDs.put(MultiSet<int>());
  SUNDANCE_VERB_HIGH(tabs << "ProdExpr: active funcs for left are: "
                     << leftFuncIDs.toString());

  SUNDANCE_VERB_HIGH(tabs << "ProdExpr: active funcs for right are: "
                     << rightFuncIDs.toString());

  int maxSpatialOrder = maxOrder(multiIndices);
  int maxDiffOrder = context.topLevelDiffOrder() + maxSpatialOrder;



  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: getting left operand's nonzeros");
  leftEvaluatable()->findNonzeros(context, miLeft,
                                  leftFuncIDs,
                                  regardFuncsAsConstant);

  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: getting right operand's nonzeros");
  rightEvaluatable()->findNonzeros(context, miRight,
                                   rightFuncIDs,
                                   regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> leftSparsity 
    = leftEvaluatable()->sparsitySubset(context, miLeft, leftFuncIDs, true);

  RefCountPtr<SparsitySubset> rightSparsity 
    = rightEvaluatable()->sparsitySubset(context, miRight, rightFuncIDs, true);



  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: left sparsity subset is " 
                       << endl << *leftSparsity);
  SUNDANCE_VERB_MEDIUM(tabs << "ProdExpr: right sparsity subset is " 
                       << endl << *rightSparsity);

  for (int i=0; i<leftSparsity->numDerivs(); i++)
    {
      const MultipleDeriv& dLeft = leftSparsity->deriv(i);
     
      for (int j=0; j<rightSparsity->numDerivs(); j++)
        {
          const MultipleDeriv& dRight = rightSparsity->deriv(j);
        
          SUNDANCE_VERB_HIGH(tabs << "ProductExpr " + toString() 
                             << " considering L=" << dLeft.toString() 
                             << " R=" << dRight.toString());

          /* Skip combinations of functional derivatives that contribute
           * only to derivatives of an order we don't need */
          if (dRight.order() + dLeft.order() > maxDiffOrder) 
            {
              SUNDANCE_VERB_HIGH(tabs << "rejected because order is higher "
                                 "than needed here");
              continue;
            }

          /* Skip combinations of spatial derivs of greater order
           * than the max order of our multiindices */
          if (dRight.spatialOrder() + dLeft.spatialOrder() > maxSpatialOrder) 
            {
              SUNDANCE_VERB_HIGH(tabs << "rejected because spatial "
                                 "order is higher than needed here");
              continue;
            }

          MultiIndex netDeriv = dRight.spatialDeriv() + dLeft.spatialDeriv();
          if (dRight.spatialOrder()==dRight.order()
              && dLeft.spatialOrder()==dLeft.order()
              && !multiIndices.contains(netDeriv))
            {
              SUNDANCE_VERB_HIGH(tabs << "rejected because spatial "
                                 "deriv combination is not needed");
              continue;
            }

              
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
          
          SUNDANCE_VERB_HIGH(tabs << "ProductExpr " + toString() 
                             << "accepting L=" << dLeft.toString() 
                             << " R=" << dRight.toString());
          SUNDANCE_VERB_HIGH(tabs << "ProductExpr " + toString() 
                             << ":  " 
                             << funcs);

          
          if (!(activeFuncIDs.contains(funcs) 
                || (funcs.size()==0 && multiIndices.contains(netDeriv))))
            {
              SUNDANCE_VERB_HIGH(tabs << "skipping " << funcs);
              continue;
            }


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


  SUNDANCE_VERB_HIGH(tabs << "ProductExpr " + toString() << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "ProductExpr " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}
