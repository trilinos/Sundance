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

#include "SundanceChainRuleEvaluator.hpp"
#include "SundanceCombinatorialUtils.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

ChainRuleEvaluator::ChainRuleEvaluator(const ExprWithChildren* expr, 
                                       const EvalContext& context)
  : SubtypeEvaluator<ExprWithChildren>(expr, context), 
    expansions_(),
    childEvaluators_(expr->numChildren()),
    childSparsity_(expr->numChildren()),
    constArgDerivMap_(),
    varArgDerivMap_(),
    zerothDerivResultIndex_(-1),
    zerothDerivIsConstant_(false)
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "ChainRuleEvaluator base class ctor for " 
                     << expr->toString());
  for (int i=0; i<numChildren(); i++)
    {
      childEvaluators_[i] = expr->evaluatableChild(i)->evaluator(context);
      childEvaluators_[i]->addClient();
      childSparsity_[i] = expr->evaluatableChild(i)->sparsitySuperset(context);
    }
}

SundanceUtils::Map<OrderedPair<int, int>, Array<Array<int> > >& ChainRuleEvaluator::compMap()
{
  static Map<OrderedPair<int, int>, Array<Array<int> > > rtn;
  return rtn;
}

void ChainRuleEvaluator::resetNumCalls() const
{
  for (int i=0; i<numChildren(); i++)
    {
      childEvaluators_[i]->resetNumCalls();
    }
  Evaluator::resetNumCalls();
}


void ChainRuleEvaluator::addConstArgDeriv(const MultiSet<int>& df, int index)
{
  constArgDerivMap_.put(df, index);
}

void ChainRuleEvaluator::addVarArgDeriv(const MultiSet<int>& df, int index)
{
  varArgDerivMap_.put(df, index);
}

int ChainRuleEvaluator::constArgDerivIndex(const MultiSet<int>& df) const
{
  TEST_FOR_EXCEPTION(!constArgDerivMap_.containsKey(df), InternalError,
                     "argument derivative " << df << " not found in constant "
                     "argument derivative map");

  return constArgDerivMap_.get(df);
}

int ChainRuleEvaluator::varArgDerivIndex(const MultiSet<int>& df) const
{
  TEST_FOR_EXCEPTION(!varArgDerivMap_.containsKey(df), InternalError,
                     "argument derivative " << df << " not found in variable "
                     "argument derivative map");

  return varArgDerivMap_.get(df);
}


const Array<Array<int> >& ChainRuleEvaluator::nComps(int N, int n) const
{
  OrderedPair<int,int> key(n,N);
  if (!compMap().containsKey(key))
    {
      compMap().put(key, compositions(N)[n-1]);
    }
  return compMap().get(key);
}


double ChainRuleEvaluator::fact(int n) const
{
  TEST_FOR_EXCEPTION(n<0, InternalError, "negative argument " << n << " to factorial");
  if (n==0 || n==1) return 1.0;
  return n*fact(n-1);
}

double ChainRuleEvaluator::choose(int N, int n) const
{
  return fact(N)/fact(n)/fact(N-n);
}

double ChainRuleEvaluator::stirling2(int n, int k) const
{
  if (n < k) return 0;
  if (n == k) return 1;
  if (k<=0) return 0;
  if (k==1) return 1;
  if (n-1 == k) return choose(n, 2);
  return k*stirling2(n-1, k) + stirling2(n-1, k-1);
}


MultipleDeriv ChainRuleEvaluator::makeMD(const Array<Deriv>& d) 
{
  MultipleDeriv rtn;
  for (unsigned int i=0; i<d.size(); i++)
    {
      rtn.put(d[i]);
    }
  return rtn;
}


Set<MultiSet<MultipleDeriv> > ChainRuleEvaluator::chainRuleBins(const MultipleDeriv& d,
                                                                const MultiSet<int>& q)
{
  int n = q.size();
  Array<Array<Array<Deriv> > > bins = binnings(d, n);

  Set<MultiSet<MultipleDeriv> > rtn;

  for (unsigned int i=0; i<bins.size(); i++)
    {
      MultiSet<MultipleDeriv> b;
      for (unsigned int j=0; j<bins[i].size(); j++)
        {
          b.put(makeMD(bins[i][j]));
        }
      rtn.put(b);
    }


  return rtn;
}


int ChainRuleEvaluator::derivComboMultiplicity(const MultiSet<MultipleDeriv>& b) const
{
  /* ugly brute force code -- there has got to be a better way to compute multiplicities */

  MultipleDeriv dTot;
  Array<MultiSet<Deriv> > derivSets(b.size());
  Array<Array<Deriv> > derivArrays(b.size());
  Set<Deriv> allDerivs;
  int k=0;
  bool allDerivsAreDistinct = true;
  bool allDerivsAreIdentical = true;
  for (MultiSet<MultipleDeriv>::const_iterator i=b.begin(); i!=b.end(); i++, k++)
    {
      for (MultipleDeriv::const_iterator j=i->begin(); j!=i->end(); j++)
        {
          derivSets[k].put(*j);
          derivArrays[k].append(*j);
          dTot.put(*j);
          if (allDerivs.contains(*j)) allDerivsAreDistinct = false;
          if (allDerivs.size()>0 && !allDerivs.contains(*j)) allDerivsAreIdentical = false;
          allDerivs.put(*j);
        }
    }
  int totOrder = dTot.order();

  /* eliminate 4th order or higher */
  TEST_FOR_EXCEPTION(totOrder > 3, InternalError,
                     "deriv order " << totOrder << " not supported");

  if (b.size()==1) return 1;  /* handles case with a single multiple deriv */
  if (totOrder == b.size()) return 1; /* handles case with N first derivs */

  /* The only remaining cases are total order = 3, grouped in 2 bins. 
   * Order=3 in 1 bin and order=3 in 3 bins have been delat with already. 
   * The ways of grouping 3 derivatives in 2 bins are:
   * {(a,bc), (b,ac), (c, ab)}.
   * If all three first derivs are the same, we have (a,aa) with multiplicity 3.
   * If all are distinct, we have (a,bc) with multiplicity 1. 
   * If two coincide, we have either
   * (a,ba) with multiplicity 2, or
   * (b,aa) with multiplicity 1.
   */
  TEST_FOR_EXCEPTION(derivArrays.size() != 2, InternalError,
                     "unexpected size=" << derivArrays.size());

  if (allDerivsAreIdentical) return 3;
  if (allDerivsAreDistinct) return 1;

  if (derivArrays[0].size()==1) 
    {
      if (derivSets[1].contains(derivArrays[0][0])) return 2;
      return 1;
    }
  else
    {
      if (derivSets[0].contains(derivArrays[1][0])) return 2;
      return 1;
    }
}


void ChainRuleEvaluator::init(const ExprWithChildren* expr, 
                              const EvalContext& context)
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "ChainRuleEvaluator::init() for " 
                    << expr->toString());

  const Set<MultipleDeriv>& C = expr->findC(context);
  const Set<MultipleDeriv>& R = expr->findR(context);

  Array<Set<MultipleDeriv> > argV(expr->numChildren());
  Array<Set<MultipleDeriv> > argC(expr->numChildren());
  Array<Set<MultipleDeriv> > argR(expr->numChildren());

  for (int i=0; i<numChildren(); i++)
    {
      argV[i] = expr->evaluatableChild(i)->findV(context);
      argC[i] = expr->evaluatableChild(i)->findC(context);
      argR[i] = expr->evaluatableChild(i)->findR(context);
    }
  SUNDANCE_VERB_HIGH(tabs << "sparsity = " << *(this->sparsity()));
  typedef Set<MultipleDeriv>::const_iterator iter;

  int count=0;
  int vecResultIndex = 0;
  int constResultIndex = 0;
  for (iter md=R.begin(); md!=R.end(); md++, count++)
    {
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "working out evaluator for " << *md);
      int N = md->order();
      bool resultIsConstant = C.contains(*md);
      int resultIndex;
      if (resultIsConstant)
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "result is constant, const index=" << constResultIndex);
          addConstantIndex(count, constResultIndex);
          resultIndex = constResultIndex;
          constResultIndex++;
        }
      else
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "result is variable, vec index=" << vecResultIndex);
          addVectorIndex(count, vecResultIndex);
          resultIndex = vecResultIndex;
          vecResultIndex++;
        }

      SUNDANCE_VERB_HIGH(tab1 << "order=" << N);
      
      if (N==0)
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "zeroth deriv index=" << resultIndex);
          zerothDerivIsConstant_ = resultIsConstant;
          zerothDerivResultIndex_ = resultIndex;
          continue;
        }


      
      RefCountPtr<ChainRuleSum> sum 
        = rcp(new ChainRuleSum(resultIndex, resultIsConstant));


      for (int n=1; n<=N; n++)
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "n=" << n);
          Array<Array<Array<Deriv> > > bins = binnings(*md, n);
          SUNDANCE_VERB_HIGH(tab2 << "bins=" << bins);
          const Set<MultiSet<int> >& QW = expr->findQ_W(n, context);
          const Set<MultiSet<int> >& QC = expr->findQ_C(n, context);
          SUNDANCE_VERB_HIGH(tab2 << "Q_W=" << QW);
          SUNDANCE_VERB_HIGH(tab2 << "Q_C=" << QC);
          for (Set<MultiSet<int> >::const_iterator 
                 j=QW.begin(); j!=QW.end(); j++)
            {
              Tabs tab3;
              SUNDANCE_VERB_HIGH(tab3 << "arg index set =" << *j);
              bool argDerivIsConstant = QC.contains(*j);
              int argDerivIndex = -1;
              if (argDerivIsConstant) argDerivIndex = constArgDerivIndex(*j);
              else argDerivIndex = varArgDerivIndex(*j);
              Array<int> J = j->elements();
              Array<DerivProduct> pSum;

              Set<MultiSet<MultipleDeriv> > bins = chainRuleBins(*md, *j);
              SUNDANCE_VERB_HIGH(tab2 << "CR bins=" << bins);

              for (Set<MultiSet<MultipleDeriv> >::const_iterator 
                     k=bins.begin(); k!=bins.end(); k++)
                {
                  Tabs tab4;
                  SUNDANCE_VERB_HIGH(tab4 << "bin set =" << *k);
                  bool termIsZero = false;
                  double multiplicity = derivComboMultiplicity(*k);
                  SUNDANCE_VERB_HIGH(tab4 << "multiplicity =" << multiplicity);
                  DerivProduct prod(multiplicity);
                  int count=0;
                  for (MultiSet<MultipleDeriv>::const_iterator 
                         i=k->begin(); i!=k->end(); i++, count++)
                    {
                      int argIndex = J[count];
                      MultipleDeriv mdArg = *i;
                      const RefCountPtr<SparsitySuperset>& argSp 
                        = childSparsity_[argIndex];
                      const RefCountPtr<Evaluator>& argEv
                        = childEvaluators_[argIndex];
                      int rawValIndex = argSp->getIndex(mdArg);
                      if (argV[argIndex].contains(mdArg))
                        {
                          int valIndex 
                            = argEv->vectorIndexMap().get(rawValIndex);
                          prod.addVariableFactor(IndexPair(argIndex, valIndex));
                        }
                      else if (argC[argIndex].contains(mdArg))
                        {
                          int valIndex 
                            = argEv->constantIndexMap().get(rawValIndex);
                          prod.addConstantFactor(IndexPair(argIndex, valIndex));
                        }
                      else
                        {
                          termIsZero = true;
                          break;
                        }
                    }
                  if (termIsZero) 
                    {
                      SUNDANCE_VERB_HIGH(tab4 << "term is zero");
                      continue;
                    }
                  SUNDANCE_VERB_HIGH(tab4 << "term is nonzero");
                  pSum.append(prod);
                }
              sum->addTerm(argDerivIndex,
                           argDerivIsConstant, 
                           pSum);
            }
        }
      TEST_FOR_EXCEPTION(sum->numTerms()==0, InternalError,
                         "Empty sum in chain rule expansion for derivative "
                         << *md);
      expansions_.append(sum);
    }

  SUNDANCE_VERB_HIGH(tabs << "num constant results: " 
                     << this->sparsity()->numConstantDerivs());

  SUNDANCE_VERB_HIGH(tabs << "num var results: " 
                     << this->sparsity()->numVectorDerivs());

  
}



void ChainRuleEvaluator::internalEval(const EvalManager& mgr,
                                      Array<double>& constantResults,
                                      Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(chainRuleEvalTimer());
  Tabs tabs;

  SUNDANCE_VERB_LOW(tabs << "ChainRuleEvaluator::eval() expr=" 
                    << expr()->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *(this->sparsity()));
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());

  SUNDANCE_VERB_HIGH(tabs << "num constant results: " 
                     << this->sparsity()->numConstantDerivs());

  SUNDANCE_VERB_HIGH(tabs << "num var results: " 
                     << this->sparsity()->numVectorDerivs());

  Array<RefCountPtr<Array<double> > > constantArgResults(numChildren());
  Array<RefCountPtr<Array<RefCountPtr<EvalVector> > > > varArgResults(numChildren());

  Array<double> constantArgDerivs;
  Array<RefCountPtr<EvalVector> > varArgDerivs;

  for (int i=0; i<numChildren(); i++)
    {
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "computing results for child #"
                         << i);
                         
      constantArgResults[i] = rcp(new Array<double>());
      varArgResults[i] = rcp(new Array<RefCountPtr<EvalVector> >());
      childEvaluators_[i]->eval(mgr, *(constantArgResults[i]), 
                                *(varArgResults[i]));
    }

  evalArgDerivs(mgr, constantArgResults, varArgResults,
                constantArgDerivs, varArgDerivs);

  for (unsigned int i=0; i<expansions_.size(); i++)
    {
      Tabs tab1;
      int resultIndex = expansions_[i]->resultIndex();
      bool isConstant = expansions_[i]->resultIsConstant();
      SUNDANCE_VERB_HIGH(tab1 << "doing expansion for deriv #" << i
                         << ", result index="
                         << resultIndex);
      if (isConstant)
        {
          expansions_[i]->evalConstant(mgr, constantArgResults, constantArgDerivs,
                                       constantResults[resultIndex]);
        }
      else
        {
          expansions_[i]->evalVar(mgr, constantArgResults, varArgResults,
                                  constantArgDerivs, varArgDerivs,
                                  vectorResults[resultIndex]);
        }
    }

  if (zerothDerivResultIndex_ >= 0)
    {
      SUNDANCE_VERB_HIGH(tabs << "processing zeroth-order deriv");
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "result index = " << zerothDerivResultIndex_);
      if (zerothDerivIsConstant_)
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "zeroth-order deriv is constant");
          constantResults[zerothDerivResultIndex_] = constantArgDerivs[0];
        }
      else
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "zeroth-order deriv is variable");
          vectorResults[zerothDerivResultIndex_] = varArgDerivs[0];
        }
    }


  if (verbosity() > VerbMedium)
    {
      Tabs tab1;
      cerr << tab1 << "chain rule results " << endl;
      this->sparsity()->print(cerr, vectorResults,
                            constantResults);
    }

  SUNDANCE_VERB_LOW(tabs << "ChainRuleEvaluator::eval() done"); 
}




namespace SundanceCore
{
  namespace Internal
  {

    MultipleDeriv makeDeriv(const Expr& a)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

      TEST_FOR_EXCEPT(aPtr==0);

      Deriv d = new FunctionalDeriv(aPtr, MultiIndex());
      MultipleDeriv rtn;
      rtn.put(d);
      return rtn;
    }


    MultipleDeriv makeDeriv(const Expr& a, const Expr& b)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

      TEST_FOR_EXCEPT(aPtr==0);

      const UnknownFuncElement* bPtr
        = dynamic_cast<const UnknownFuncElement*>(b[0].ptr().get());

      TEST_FOR_EXCEPT(bPtr==0);

      Deriv da = new FunctionalDeriv(aPtr, MultiIndex());
      Deriv db = new FunctionalDeriv(bPtr, MultiIndex());
      MultipleDeriv rtn;
      rtn.put(da);
      rtn.put(db);
      return rtn;
    }



    MultipleDeriv makeDeriv(const Expr& a, const Expr& b, const Expr& c)
    {
      const UnknownFuncElement* aPtr
        = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

      TEST_FOR_EXCEPT(aPtr==0);

      const UnknownFuncElement* bPtr
        = dynamic_cast<const UnknownFuncElement*>(b[0].ptr().get());

      TEST_FOR_EXCEPT(bPtr==0);

      const UnknownFuncElement* cPtr
        = dynamic_cast<const UnknownFuncElement*>(c[0].ptr().get());

      TEST_FOR_EXCEPT(cPtr==0);

      Deriv da = new FunctionalDeriv(aPtr, MultiIndex());
      Deriv db = new FunctionalDeriv(bPtr, MultiIndex());
      Deriv dc = new FunctionalDeriv(cPtr, MultiIndex());
      MultipleDeriv rtn;
      rtn.put(da);
      rtn.put(db);
      rtn.put(dc);
      return rtn;
    }
  }
}
