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
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
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
                         << *(argSparsitySubset()));

    Map<const DiscreteFuncElementEvaluator*, int> funcToIndexMap;

    int vecResultIndex = 0;
    int constResultIndex = 0;
  
    for (int i=0; i<this->sparsity()->numDerivs(); i++)
      {
        if (this->sparsity()->state(i)==ConstantDeriv)
          {
            addConstantIndex(i, constResultIndex);
            resultIndices_[i] = constResultIndex++;
            isConstant_[i] = true;
          }
        else
          {
            addVectorIndex(i, vecResultIndex);
            resultIndices_[i] = vecResultIndex++;
            isConstant_[i] = false;
          }
      }


  
    SUNDANCE_VERB_MEDIUM(tab0 << "looking for my derivs induced by arg derivs");
    for (int i=0; i<argSparsitySubset()->numDerivs(); i++)
      {
        Tabs tab1;
        Map<MultipleDeriv, DerivState> monomials;
        Map<MultipleDeriv, Deriv> funcTerms;
        expr->getResultDerivs(argSparsitySubset()->deriv(i), 
                              argSparsitySubset()->state(i),
                              monomials, 
                              funcTerms);
      
        SUNDANCE_VERB_MEDIUM(tab1 << "looking for effect of argument deriv " 
                             << argSparsitySubset()->deriv(i));

        
        {
          Tabs tab2;
          SUNDANCE_VERB_MEDIUM(tab2 << "found monomials " << monomials
                               << " and function terms " << funcTerms);
      
          /* we'll need the index of the argument's deriv in the superset
           * of results */
          int argIndexInSuperset = argSparsitySuperset()->getIndex(argSparsitySubset()->deriv(i));
          TEST_FOR_EXCEPTION(argIndexInSuperset==-1, InternalError,
                             "derivative " << argSparsitySubset()->deriv(i)
                             << " found in arg subset but not in arg superset");

          SUNDANCE_VERB_MEDIUM(tab2 << "looking for effects of function terms");

          for (Map<MultipleDeriv, Deriv>::const_iterator 
                 iter = funcTerms.begin(); iter != funcTerms.end(); iter++)
            {
              Tabs tab3;
              const MultipleDeriv& target = iter->first;
              const Deriv& d = iter->second;

              SUNDANCE_VERB_MEDIUM(tab3 << "derivative " << d << " contributes to "
                                   "functional deriv " << target);
          

              TEST_FOR_EXCEPTION(d.isCoordDeriv(), InternalError,
                                 "Coord deriv detected where a functional "
                                 "deriv was expected in "
                                 "DiffOpEvaluator ctor");
              const FuncElementBase* f = d.funcDeriv()->func();
              const MultiIndex& mi = d.funcDeriv()->multiIndex();
          
          
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
          
              bool coeffIsConstant = argSparsitySubset()->state(i)==ConstantDeriv;
              
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

              int resultIndex = this->sparsity()->getIndex(target);

              if (coeffIsConstant)
                {
                  Tabs tab4;
                  SUNDANCE_VERB_MEDIUM(tab4 << "found constant coeff function term: arg deriv "
                                       << argSparsitySubset()->deriv(i) 
                                       << " times "
                                       << d);
                  int constArgIndex = argEval()->constantIndexMap().get(argIndexInSuperset);
                  constantCoeffFuncIndices_[resultIndex].append(fIndex);
                  constantCoeffFuncMi_[resultIndex].append(miIndex);
                  constantFuncCoeffs_[resultIndex].append(constArgIndex);
                }
              else
                {     
                  Tabs tab4;
                  SUNDANCE_VERB_MEDIUM(tab4 << "found vec coeff function term: arg deriv "
                                       << argSparsitySubset()->deriv(i) 
                                       << " times "
                                       << d);
                  int vectorArgIndex = argEval()->vectorIndexMap().get(argIndexInSuperset);
                  vectorCoeffFuncIndices_[resultIndex].append(fIndex);
                  vectorCoeffFuncMi_[resultIndex].append(miIndex);
                  vectorFuncCoeffs_[resultIndex].append(vectorArgIndex);
                }

            }
    
      
          SUNDANCE_VERB_MEDIUM(tab1 << "looking for effects of monomials");
      
          for (Map<MultipleDeriv, DerivState>::const_iterator 
                 iter = monomials.begin(); iter != monomials.end(); 
               iter++)
            {
              Tabs tab3;
              int resultIndex = this->sparsity()->getIndex(iter->first);
              SUNDANCE_VERB_MEDIUM(tab3 << "monomial " << iter->first
                                   << " contributes to result index " 
                                   << resultIndex);
              if (iter->second==ConstantDeriv)
                {
                  Tabs tab4;
                  SUNDANCE_VERB_MEDIUM(tab4 << "found constant-valued monomial: "
                                       "arg deriv " 
                                       << argSparsitySubset()->deriv(i));
                  int constArgIndex = argEval()->constantIndexMap().get(argIndexInSuperset);
                  constantMonomials_[resultIndex].append(constArgIndex);
                }
              else
                {
                  Tabs tab4;
                  SUNDANCE_VERB_MEDIUM(tab4 << "found vector-valued monomial: "
                                       "arg deriv " << argSparsitySubset()->deriv(i));
                  int vectorArgIndex = argEval()->vectorIndexMap().get(argIndexInSuperset);
                  vectorMonomials_[resultIndex].append(vectorArgIndex);
                }
            }
        }
      }
  }
  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "instruction tables for summing chain rule" << endl;
      for (int i=0; i<this->sparsity()->numDerivs(); i++)
        {
          Tabs tab1;
          cerr << tab1 << "deriv " << i << endl;
          {
            Tabs tab2;
            cerr << tab2 << "constant monomials: " << constantMonomials_[i]
                 << endl;
            cerr << tab2 << "vector monomials: " << vectorMonomials_[i]
                 << endl;
            
            cerr << tab2 << "constant coeff functions: " << endl;
            for (unsigned int j=0; j<constantFuncCoeffs_[i].size(); j++)
              {
                Tabs tab3;
                cerr << tab3 << "func=" << constantCoeffFuncIndices_[i][j]
                     << " mi=" << constantCoeffFuncMi_[i][j] << endl;
              } 
            cerr << tab2 << "vector coeff functions: " << endl;
            for (unsigned int j=0; j<vectorFuncCoeffs_[i].size(); j++)
              {
                Tabs tab3;
                cerr << tab3 << "func=" << vectorCoeffFuncIndices_[i][j]
                     << " mi=" << vectorCoeffFuncMi_[i][j] << endl;
              }
            
          }
        }
    }
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
      cerr << tabs << "DiffOp operand results" << endl;
      argSparsitySuperset()->print(cerr, argVectorResults,
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
      SUNDANCE_VERB_MEDIUM(tab1 << "adding vector monomials");
      const Array<int>& vm = vectorMonomials_[i];
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
      SUNDANCE_VERB_MEDIUM(tab1 << "adding func terms with constant coeffs");
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
      SUNDANCE_VERB_MEDIUM(tab1 << "adding func terms with vector coeffs");
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
      cerr << tabs << "operand results" << endl;
      sparsity()->print(cerr, vectorResults,
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
