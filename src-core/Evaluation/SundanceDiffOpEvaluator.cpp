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

#include "SundanceDiffOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceZeroExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;




DiffOpEvaluator
::DiffOpEvaluator(const DiffOp* expr,
                  const EvalContext& context)
  : UnaryEvaluator<DiffOp>(expr, context),
    isConstant_(this->sparsity()->numDerivs()),
    resultIndices_(this->sparsity()->numDerivs()),
    constantMonomials_(this->sparsity()->numDerivs()),
    vectorMonomials_(this->sparsity()->numDerivs()),
    constantFuncCoeffs_(this->sparsity()->numDerivs()),
    vectorFuncCoeffs_(this->sparsity()->numDerivs()),
    funcEvaluators_(),
    constantCoeffFuncIndices_(this->sparsity()->numDerivs()),
    constantCoeffFuncMi_(this->sparsity()->numDerivs()),
    vectorCoeffFuncIndices_(this->sparsity()->numDerivs()),
    vectorCoeffFuncMi_(this->sparsity()->numDerivs())
{
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing diff op evaluator for " 
                    << expr->toString());

  {
    Tabs tab0;
  
    SUNDANCE_VERB_MEDIUM(tab0 << "return sparsity " << endl << *(this->sparsity)());

    SUNDANCE_VERB_MEDIUM(tab0 << "argument sparsity subset " << endl 
                         << *(argSparsitySuperset()));

    Map<const DiscreteFuncElementEvaluator*, int> funcToIndexMap;

    int vecResultIndex = 0;
    int constResultIndex = 0;
  
    for (int i=0; i<this->sparsity()->numDerivs(); i++)
      {
        Tabs tab1;
        const MultipleDeriv& resultDeriv = this->sparsity()->deriv(i);
        SUNDANCE_VERB_HIGH(tab0 << "working out procedure for computing " 
                           << resultDeriv);

        if (this->sparsity()->state(i)==ConstantDeriv)
          {
            Tabs tab2;
            addConstantIndex(i, constResultIndex);
            resultIndices_[i] = constResultIndex++;
            isConstant_[i] = true;
            SUNDANCE_VERB_HIGH(tab2 << "deriv is constant, will be stored at index "
                               << resultIndices_[i] << " in the const result array");
          }
        else
          {
            Tabs tab2;
            addVectorIndex(i, vecResultIndex);
            resultIndices_[i] = vecResultIndex++;
            isConstant_[i] = false;
            SUNDANCE_VERB_HIGH(tab2 << "deriv is variable, will be stored at index "
                               << resultIndices_[i] << " in the var result array");
          }

        int order = resultDeriv.order();
        const Set<MultipleDeriv>& RArg 
          = argExpr()->findR(order, context);
        const Set<MultipleDeriv>& RArgPlus
          = argExpr()->findR(order+1, context);
        const Set<MultipleDeriv>& W1Arg 
          = argExpr()->findW(1, context);

        
        SUNDANCE_VERB_HIGH(tab1 << "RArg = " << RArg);
        SUNDANCE_VERB_HIGH(tab1 << "RArgPlus = " << RArgPlus);
        SUNDANCE_VERB_HIGH(tab1 << "W1Arg = " << W1Arg);

        Set<MultipleDeriv> funcTermCoeffs 
          = RArgPlus.intersection(increasedDerivs(resultDeriv, W1Arg));
        SUNDANCE_VERB_HIGH(tab1 << "function term coeffs = " << funcTermCoeffs);

        SUNDANCE_VERB_HIGH(tab1 << "getting direct chain rule terms");


        for (Set<MultipleDeriv>::const_iterator 
               j=funcTermCoeffs.begin(); j != funcTermCoeffs.end(); j++)
          {
            Tabs tab2;
            SUNDANCE_VERB_HIGH(tab2 << "getting coefficient of " << *j);

            int argIndex = argSparsitySuperset()->getIndex(*j);
            TEST_FOR_EXCEPTION(argIndex==-1, RuntimeError,
                               "Derivative " << *j << " expected in argument "
                               "but not found");

            Deriv lambda = remainder(*j, resultDeriv);

            if (lambda.isCoordDeriv())
              {
                Tabs tab3;
                SUNDANCE_VERB_HIGH(tab2 << "detected coordinate deriv");
                if (lambda.coordDeriv()->dir()!=expr->mi().firstOrderDirection())
                  {
                    SUNDANCE_VERB_HIGH(tab2 << "direction mismatch, skipping");
                    continue;
                  }
                const DerivState& argState = argSparsitySuperset()->state(argIndex);
                if (argState==ConstantDeriv)
                  {
                    int constArgIndex = argEval()->constantIndexMap().get(argIndex);
                    constantMonomials_[i].append(constArgIndex);
                  }
                else
                  {
                    int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
                    vectorMonomials_[i].append(vectorArgIndex);
                  }
              }
            else
              {
                Tabs tab3;
                SUNDANCE_VERB_HIGH(tab3 << "detected functional deriv " << lambda);
                const FuncElementBase* f = lambda.funcDeriv()->func();
                const MultiIndex& mi = expr->mi() + lambda.funcDeriv()->multiIndex(); 
                SUNDANCE_VERB_HIGH(tab3 << "modified multiIndex is " << mi);

                const TestFuncElement* t 
                  = dynamic_cast<const TestFuncElement*>(f);
                if (t != 0) continue;

                const UnknownFuncElement* u 
                  = dynamic_cast<const UnknownFuncElement*>(f);
                TEST_FOR_EXCEPTION(u==0, InternalError,
                                   "Non-unknown function detected where an unknown "
                                   "function was expected in "
                                   "DiffOpEvaluator ctor");


                const EvaluatableExpr* evalPt = u->evalPt();
                const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(evalPt);
                if (z != 0) continue;
                TEST_FOR_EXCEPTION(z != 0, InternalError,
                                   "DiffOpEvaluator detected identically zero "
                                   "function");

                const DiscreteFuncElement* df 
                  = dynamic_cast<const DiscreteFuncElement*>(evalPt);
          
                TEST_FOR_EXCEPTION(df==0, InternalError,
                                   "DiffOpEvaluator ctor: evaluation point of "
                                   "unknown function " << u->toString() 
                                   << " is not a discrete function");

                const SymbolicFuncElementEvaluator* uEval 
                  = dynamic_cast<const SymbolicFuncElementEvaluator*>(u->evaluator(context).get());

                const DiscreteFuncElementEvaluator* dfEval = uEval->dfEval();


                TEST_FOR_EXCEPTION(dfEval==0, InternalError,
                                   "DiffOpEvaluator ctor: evaluator for "
                                   "evaluation point is not a "
                                   "DiscreteFuncElementEvaluator");

                TEST_FOR_EXCEPTION(!dfEval->hasMultiIndex(mi), InternalError,
                                   "DiffOpEvaluator ctor: evaluator for "
                                   "discrete function " << df->toString()
                                   << " does not know about multiindex "
                                   << mi.toString());
          
                int fIndex;
                int miIndex = dfEval->miIndex(mi);
          
                if (funcToIndexMap.containsKey(dfEval))
                  {
                    fIndex = funcToIndexMap.get(dfEval);
                  }
                else
                  {
                    fIndex = funcEvaluators_.size();
                    funcEvaluators_.append(dfEval);
                    funcToIndexMap.put(dfEval, fIndex);
                  }

            
                const DerivState& argState = argSparsitySuperset()->state(argIndex);
                if (argState==ConstantDeriv)
                  {
                    int constArgIndex = argEval()->constantIndexMap().get(argIndex);
                    constantCoeffFuncIndices_[i].append(fIndex);
                    constantCoeffFuncMi_[i].append(miIndex);
                    constantFuncCoeffs_[i].append(constArgIndex);
                  }
                else
                  {
                    int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
                    vectorCoeffFuncIndices_[i].append(fIndex);
                    vectorCoeffFuncMi_[i].append(miIndex);
                    vectorFuncCoeffs_[i].append(vectorArgIndex);
                  }
              }
          }


        SUNDANCE_VERB_HIGH(tab1 << "getting indirect chain rule terms");
        Set<MultipleDeriv> isolatedTerms 
          = RArg.intersection(backedDerivs(resultDeriv, W1Arg));
        SUNDANCE_VERB_HIGH(tab1 << "isolated terms = " << isolatedTerms);

        for (Set<MultipleDeriv>::const_iterator 
               j=isolatedTerms.begin(); j != isolatedTerms.end(); j++)
          {
            int argIndex = argSparsitySuperset()->getIndex(*j);
            TEST_FOR_EXCEPTION(argIndex==-1, RuntimeError,
                               "Derivative " << *j << " expected in argument "
                               "but not found");
            const DerivState& argState = argSparsitySuperset()->state(argIndex);
            if (argState==ConstantDeriv)
              {
                int constArgIndex = argEval()->constantIndexMap().get(argIndex);
                constantMonomials_[i].append(constArgIndex);
              }
            else
              {
                int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
                vectorMonomials_[i].append(vectorArgIndex);
              }
          }
        /*
        TEST_FOR_EXCEPTION(
          constantMonomials_[i].size()==0U
          && vectorMonomials_[i].size()==0U
          && constantFuncCoeffs_[i].size()==0U
          && vectorFuncCoeffs_[i].size()==0U,
          RuntimeError,
          "no instructions found for deriv " << resultDeriv);
        */
      }
  }

  if (verbosity() > VerbMedium)
    {
      Out::os() << tabs << "instruction tables for summing chain rule" << endl;
      for (int i=0; i<this->sparsity()->numDerivs(); i++)
        {
          Tabs tab1;
          Out::os() << tab1 << "deriv " << sparsity()->deriv(i) << endl;
          {
            Tabs tab2;
            Out::os() << tab2 << "constant monomials: " << constantMonomials_[i]
                 << endl;
            Out::os() << tab2 << "vector monomials: " << vectorMonomials_[i]
                 << endl;
            
            Out::os() << tab2 << "constant coeff functions: " << endl;
            for (unsigned int j=0; j<constantFuncCoeffs_[i].size(); j++)
              {
                Tabs tab3;
                Out::os() << tab3 << "func=" << constantCoeffFuncIndices_[i][j]
                     << " mi=" << constantCoeffFuncMi_[i][j] << endl;
              } 
            Out::os() << tab2 << "vector coeff functions: " << endl;
            for (unsigned int j=0; j<vectorFuncCoeffs_[i].size(); j++)
              {
                Tabs tab3;
                Out::os() << tab3 << "func=" << vectorCoeffFuncIndices_[i][j]
                     << " mi=" << vectorCoeffFuncMi_[i][j] << endl;
              }
            
          }
        }
    }
}


Deriv DiffOpEvaluator::remainder(const MultipleDeriv& big, 
                                 const MultipleDeriv& little) const 
{
  Tabs tab;
  SUNDANCE_VERB_EXTREME(tab << "computing remainder: big=" << big << ", little="
                        << little);
  TEST_FOR_EXCEPT(big.order()-little.order() != 1);

  MultipleDeriv r;
  if (little.order()==0) r = big;
  else r = big.factorOutDeriv(little);

  SUNDANCE_VERB_EXTREME(tab << "remainder = " << r);

  TEST_FOR_EXCEPT(r.order() != 1);

  return *(r.begin());
}

Set<MultipleDeriv> DiffOpEvaluator
::increasedDerivs(const MultipleDeriv& mu,
                  const Set<MultipleDeriv>& W1) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "computing increased derivs");
  Set<MultipleDeriv> rtn;
  for (Set<MultipleDeriv>::const_iterator i=W1.begin(); i!=W1.end(); i++)
    {
      MultipleDeriv md = *i;
      TEST_FOR_EXCEPT(md.order() != 1);
      Deriv lambda = *(md.begin());
      MultipleDeriv lambdaMu = mu;
      lambdaMu.put(lambda);
      rtn.put(lambdaMu);
    }
  SUNDANCE_VERB_HIGH(tabs << "increased derivs = " << rtn);
  return rtn;
}

Set<MultipleDeriv> DiffOpEvaluator
::backedDerivs(const MultipleDeriv& mu,
               const Set<MultipleDeriv>& W1) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "computing backed-out derivs for mu= " << mu
                     << ", W1=" << W1);
  Set<MultipleDeriv> rtn;
  if (mu.order() != 0) 
    {
      const MultiIndex& alpha = expr()->mi();

      for (Set<MultipleDeriv>::const_iterator i=W1.begin(); i!=W1.end(); i++)
        {
          const MultipleDeriv& md = *i;
          TEST_FOR_EXCEPT(md.order() != 1);
          Deriv lambda = *(md.begin());
          if (lambda.isCoordDeriv()) continue;
          int lambda_fid = lambda.funcDeriv()->funcComponentID();
          const MultiIndex& lambda_mi = lambda.funcDeriv()->multiIndex(); 
          for (MultipleDeriv::const_iterator j=mu.begin(); j!=mu.end(); j++)
            {
              const Deriv& d = *j;
              if (d.isCoordDeriv()) continue;
              int d_fid = d.funcDeriv()->funcComponentID();
              const MultiIndex& d_mi = d.funcDeriv()->multiIndex(); 
              if (d_fid != lambda_fid) continue;
              if (!(alpha + lambda_mi == d_mi)) continue;
              MultipleDeriv z = mu.factorOutDeriv(d);
              z.put(lambda);
              rtn.put(z);
            }
        }
    }
  SUNDANCE_VERB_HIGH(tabs << "backed-out derivs = " << rtn);
  return rtn;
}



void DiffOpEvaluator::internalEval(const EvalManager& mgr,
                                   Array<double>& constantResults,
                                   Array<RefCountPtr<EvalVector> >& vectorResults)  const
{
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "DiffOpEvaluator::eval() expr=" 
                    << expr()->toString());

  /* evaluate the argument */
  Array<RefCountPtr<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  evalOperand(mgr, argConstantResults, argVectorResults);


  if (verbosity() > VerbLow)
    {
      Out::os() << tabs << "DiffOp operand results" << endl;
      argSparsitySuperset()->print(Out::os(), argVectorResults,
                                   argConstantResults);
    }



  /* evaluate the required discrete functions */
  SUNDANCE_VERB_MEDIUM(tabs << "evaluating discrete functions");

  Array<Array<RefCountPtr<EvalVector> > > funcVectorResults(funcEvaluators_.size());
  Array<double> funcConstantResults;
  for (unsigned int i=0; i<funcEvaluators_.size(); i++)
    {
      funcEvaluators_[i]->eval(mgr, funcConstantResults, funcVectorResults[i]);
    }
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());
  
  SUNDANCE_VERB_MEDIUM(tabs << "summing spatial/functional chain rule");

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      Tabs tab1;
      SUNDANCE_VERB_MEDIUM(tab1 << "working on deriv " 
                           << this->sparsity()->deriv(i));

      /* add constant monomials */
      SUNDANCE_VERB_MEDIUM(tab1 << "have " <<  constantMonomials_[i].size()
        << " constant monomials");
      double constantVal = 0.0;
      for (unsigned int j=0; j<constantMonomials_[i].size(); j++)
        {
          SUNDANCE_VERB_MEDIUM(tab1 << "adding in constant monomial (index "
                               << constantMonomials_[i][j] 
                               << " in arg results)");
          constantVal += argConstantResults[constantMonomials_[i][j]];
        }
      if (isConstant_[i])
        {
          constantResults[resultIndices_[i]] = constantVal;
          SUNDANCE_VERB_MEDIUM(tab1 << "result is constant: value=" 
                               << constantVal);
          continue;
        }

      RefCountPtr<EvalVector> result;
      bool vecHasBeenAllocated = false;

      /* add in the vector monomials */
      const Array<int>& vm = vectorMonomials_[i];
      SUNDANCE_VERB_MEDIUM(tab1 << "have " << vm.size() 
        << " vector monomials");
      for (unsigned int j=0; j<vm.size(); j++)
        {
          Tabs tab2;

          const RefCountPtr<EvalVector>& v = argVectorResults[vm[j]];

          SUNDANCE_VERB_MEDIUM(tab2 << "found term " << v->str());

          /* if we've not yet allocated a vector for the results, 
           * do so now, and set it to the initial value */ 
          if (!vecHasBeenAllocated)
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "allocated result vector");
              result = mgr.popVector();
              vecHasBeenAllocated = true;
              if (isZero(constantVal))
                {
                  result->setTo_V(v.get());
                }
              else
                {
                  result->setTo_S_add_V(constantVal, v.get());
                }
            }
          else
            {
              result->add_V(v.get());
            }
          SUNDANCE_VERB_MEDIUM(tab2 << "result is " << result->str());
        }
      
      /* add in the function terms with constant coeffs */
      const Array<int>& cf = constantFuncCoeffs_[i];
      SUNDANCE_VERB_MEDIUM(tab1 << "adding " << cf.size()
        << " func terms with constant coeffs");
      for (unsigned int j=0; j<cf.size(); j++)
        {
          Tabs tab2;
          const double& coeff = argConstantResults[cf[j]];
          int fIndex = constantCoeffFuncIndices_[i][j];
          int miIndex = constantCoeffFuncMi_[i][j];
          const RefCountPtr<EvalVector>& fValue 
            = funcVectorResults[fIndex][miIndex];

          SUNDANCE_VERB_MEDIUM(tab2 << "found term: coeff= " 
                               << coeff << ", func value=" << fValue->str());
          
          /* if we've not yet allocated a vector for the results, 
           * do so now, and set it to the initial value */ 
          if (!vecHasBeenAllocated)
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "allocated result vector");
              result = mgr.popVector();
              vecHasBeenAllocated = true;
              if (isOne(coeff))
                {
                  if (isZero(constantVal))
                    {
                      result->setTo_V(fValue.get());
                    }
                  else
                    {
                      result->setTo_S_add_V(constantVal, fValue.get());
                    }
                }
              else
                {
                  if (isZero(constantVal))
                    {
                      result->setTo_SV(coeff, fValue.get());
                    }
                  else
                    {
                      result->setTo_S_add_SV(constantVal, coeff, fValue.get());
                    }
                }
            }
          else
            {
              if (isOne(coeff))
                {
                  result->add_V(fValue.get());
                }
              else
                {
                  result->add_SV(coeff, fValue.get());
                }
            }
          SUNDANCE_VERB_MEDIUM(tab2 << "result is " << result->str());
        }

      
      /* add in the function terms with vector coeffs */
      const Array<int>& vf = vectorFuncCoeffs_[i];
      SUNDANCE_VERB_MEDIUM(tab1 << "adding " << vf.size()
        << " func terms with vector coeffs");
      for (unsigned int j=0; j<vf.size(); j++)
        {
          Tabs tab2;

          const RefCountPtr<EvalVector>& coeff = argVectorResults[vf[j]];
          int fIndex = vectorCoeffFuncIndices_[i][j];
          int miIndex = vectorCoeffFuncMi_[i][j];
          const RefCountPtr<EvalVector>& fValue 
            = funcVectorResults[fIndex][miIndex];

          SUNDANCE_VERB_MEDIUM(tab2 << "found term: coeff= " 
                               << coeff->str() << ", func value=" 
                               << fValue->str());
          
          /* if we've not yet allocated a vector for the results, 
           * do so now, and set it to the initial value */ 
          if (!vecHasBeenAllocated)
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "allocated result vector");
              result = mgr.popVector();
              vecHasBeenAllocated = true;
              result->setTo_VV(coeff.get(), fValue.get());
            }
          else
            {
              result->add_VV(coeff.get(), fValue.get());
            }
          SUNDANCE_VERB_MEDIUM(tab2 << "result is " << result->str());
        }

      TEST_FOR_EXCEPTION(!vecHasBeenAllocated, InternalError,
                         "created empty vector in DiffOpEvaluator::internalEval");
      vectorResults[resultIndices_[i]] = result;
    }

  if (verbosity() > VerbLow)
    {
      Out::os() << tabs << "operand results" << endl;
      sparsity()->print(Out::os(), vectorResults,
                        constantResults);
    }
  SUNDANCE_VERB_MEDIUM(tabs << "done chain rule");
}



void DiffOpEvaluator::resetNumCalls() const 
{
  argEval()->resetNumCalls();
  for (unsigned int i=0; i<funcEvaluators_.size(); i++) 
    {
      funcEvaluators_[i]->resetNumCalls();
    }
  Evaluator::resetNumCalls();
}
