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

#include "SundanceSumEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;



SumEvaluator::SumEvaluator(const SumExpr* se,
                           const EvalContext& context)
  : BinaryEvaluator<SumExpr>(se, context), 
    sign_(se->sign()),
    singleRightConstant_(),
    singleRightVector_(),
    singleLeftConstant_(),
    singleLeftVector_(),
    ccSums_(),
    cvSums_(),
    vcSums_(),
    vvSums_()
{
  Tabs tabs;

  if (verbosity() > VerbLow)
    {
      cerr << tabs << "initializing sum evaluator for " 
           << se->toString() << endl;
    }

  int constantCounter = 0;
  int vectorCounter = 0;

  if (verbosity() > VerbMedium)
    {
      cerr << endl << tabs << "return sparsity ";
      sparsity()->print(cerr);
      cerr << endl << "left sparsity ";
      leftSparsity()->print(cerr);
      cerr << endl << "right sparsity " ;
      rightSparsity()->print(cerr);
      cerr << endl;
      
      cerr << "left vector index map " << leftEval()->vectorIndexMap() << endl;
      cerr << "right vector index map " << rightEval()->vectorIndexMap() << endl;
      
      cerr << "left constant index map " << leftEval()->constantIndexMap() << endl;
      cerr << "right constant index map " << rightEval()->constantIndexMap() << endl;
    }

  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity()->deriv(i);
      TEST_FOR_EXCEPTION(!leftSparsity()->containsDeriv(d) 
                         && !rightSparsity()->containsDeriv(d),
                         InternalError,
                         "deriv " << d.toString() 
                         << " was not found in either left or right operand "
                         "of expr " << expr()->toString());
      
      int iLeft = -1;
      int iRight = -1;
      
      if (leftSparsity()->containsDeriv(d)) 
        {
          iLeft = leftSparsity()->getIndex(d);
        }

      if (rightSparsity()->containsDeriv(d)) 
        {
          iRight = rightSparsity()->getIndex(d);
        }

      SUNDANCE_VERB_MEDIUM(tabs << "deriv " << d);

      if (iLeft == -1) /* case where the left operand is zero */
        {
          SUNDANCE_VERB_HIGH(tabs << "left operand is zero ");
          if (rightSparsity()->state(iRight)==ConstantDeriv)
            {
              int rc = rightEval()->constantIndexMap().get(iRight);
              singleRightConstant_.append(tuple(constantCounter, rc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_VERB_HIGH(tabs << "single right constant " 
                                 << singleRightConstant_[singleRightConstant_.size()-1]);
            }
          else
            {
              int rv = rightEval()->vectorIndexMap().get(iRight);
              singleRightVector_.append(tuple(vectorCounter, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_VERB_HIGH(tabs << "single right vector " 
                                 << singleRightVector_[singleRightVector_.size()-1]);
            }
        }
      else if (iRight == -1) /* case where the right operand is zero */
        {
          SUNDANCE_VERB_HIGH(tabs << "right operand is zero ");
          if (leftSparsity()->state(iLeft)==ConstantDeriv)
            {
              int lc = leftEval()->constantIndexMap().get(iLeft);
              singleLeftConstant_.append(tuple(constantCounter, lc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_VERB_HIGH(tabs << "single left constant " 
                                 << singleLeftConstant_[singleLeftConstant_.size()-1]);
            }
          else
            {
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              singleLeftVector_.append(tuple(vectorCounter, lv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_VERB_HIGH(tabs << "single left vector " 
                                 << singleLeftVector_[singleLeftVector_.size()-1]);
            }
        }
      else /* both are nonzero */
        {
          SUNDANCE_VERB_HIGH(tabs << "both operands are nonzero ");
          bool leftIsConstant = leftSparsity()->state(iLeft)==ConstantDeriv;
          bool rightIsConstant = rightSparsity()->state(iRight)==ConstantDeriv;
          
          if (leftIsConstant && rightIsConstant)
            {
              SUNDANCE_VERB_HIGH(tabs << "both operands are constant");
              int lc = leftEval()->constantIndexMap().get(iLeft);
              int rc = rightEval()->constantIndexMap().get(iRight);
              ccSums_.append(tuple(constantCounter, lc, rc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_VERB_HIGH(tabs << "c-c sum " << ccSums_[ccSums_.size()-1]);
            }
          else if (leftIsConstant)
            {
              SUNDANCE_VERB_HIGH(tabs << "left operand is constant");
              int lc = leftEval()->constantIndexMap().get(iLeft);
              int rv = rightEval()->vectorIndexMap().get(iRight);
              cvSums_.append(tuple(vectorCounter, lc, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_VERB_HIGH(tabs << "c-v sum " << cvSums_[cvSums_.size()-1]);
            }
          else if (rightIsConstant)
            {
              SUNDANCE_VERB_HIGH(tabs << "right operand is constant");
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              int rc = rightEval()->constantIndexMap().get(iRight);
              vcSums_.append(tuple(vectorCounter, lv, rc));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_VERB_HIGH(tabs << "v-c sum " << vcSums_[vcSums_.size()-1]);
            }
          else
            {
              SUNDANCE_VERB_HIGH(tabs << "both operands are vectors");
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              int rv = rightEval()->vectorIndexMap().get(iRight);
              vvSums_.append(tuple(vectorCounter, lv, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_VERB_HIGH(tabs << "v-v sum " << vvSums_[vvSums_.size()-1]);
            }
        }
    }
}

void SumEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const 
{ 
  TimeMonitor timer(evalTimer());
  Tabs tabs;

  SUNDANCE_OUT(verbosity() > VerbSilent,
               tabs << "SumEvaluator::eval() expr=" << expr()->toString());

  /* evaluate the children */
  Array<RefCountPtr<EvalVector> > leftVectorResults; 
  Array<RefCountPtr<EvalVector> > rightVectorResults; 
  Array<double> leftConstResults;
  Array<double> rightConstResults;
  evalChildren(mgr, leftConstResults, leftVectorResults,
               rightConstResults, rightVectorResults);

  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "left operand " << endl;
      leftSparsity()->print(cerr, leftVectorResults,
                            leftConstResults);
      cerr << tabs << "right operand " << endl;
      rightSparsity()->print(cerr, rightVectorResults,
                             rightConstResults);
    }
  
  constantResults.resize(sparsity()->numConstantDerivs());
  vectorResults.resize(sparsity()->numVectorDerivs());

  /* Do constant terms with left=0 */
  for (unsigned int i=0; i<singleRightConstant_.size(); i++)
    {
      Tabs tab1;
      constantResults[singleRightConstant_[i][0]]
        = sign_*rightConstResults[singleRightConstant_[i][1]];
      SUNDANCE_VERB_MEDIUM(tab1 << "sum for "
                           << constantResultDeriv(singleRightConstant_[i][0])
                           << ": L=0, R=" 
                           << sign_*rightConstResults[singleRightConstant_[i][1]]);
    }

  /* Do constant terms with right=0 */
  for (unsigned int i=0; i<singleLeftConstant_.size(); i++)
    {
      Tabs tab1;
      constantResults[singleLeftConstant_[i][0]]
        = leftConstResults[singleLeftConstant_[i][1]];
      SUNDANCE_VERB_MEDIUM(tab1 << "sum for " 
                           << constantResultDeriv(singleLeftConstant_[i][0])
                           << ": L=" 
                           << leftConstResults[singleLeftConstant_[i][1]]
                           << " R=0");
    }

  /* Do vector terms with left=0 */
  for (unsigned int i=0; i<singleRightVector_.size(); i++)
    {
      Tabs tab1;
      if (sign_ < 0.0) rightVectorResults[singleRightVector_[i][1]]->multiply_S(sign_);
      vectorResults[singleRightVector_[i][0]]
        = rightVectorResults[singleRightVector_[i][1]];
      SUNDANCE_VERB_MEDIUM(tab1 << "sum for "
                           << vectorResultDeriv(singleRightVector_[i][0])
                           << ": L=0, R=" 
                           << sign_ << "*" << 
                           rightVectorResults[singleRightVector_[i][1]]->str());
    }

  /* Do vector terms with right=0 */
  for (unsigned int i=0; i<singleLeftVector_.size(); i++)
    { 
      Tabs tab1;
      vectorResults[singleLeftVector_[i][0]]
        = leftVectorResults[singleLeftVector_[i][1]];
      SUNDANCE_VERB_MEDIUM(tab1 << "sum for " 
                           << vectorResultDeriv(singleLeftVector_[i][0])
                           << ": L=" 
                           << leftVectorResults[singleLeftVector_[i][1]]->str()
                           << " R=0");
    }

  /** Do constant-constant terms */
  for (unsigned int i=0; i<ccSums_.size(); i++)
    {
      Tabs tab1;
      constantResults[ccSums_[i][0]]
        = leftConstResults[ccSums_[i][1]] 
        + sign_*rightConstResults[ccSums_[i][2]];
      SUNDANCE_VERB_MEDIUM(tab1 << "c-c sum for " 
                           << vectorResultDeriv(ccSums_[i][0])
                           << ": L=" << leftConstResults[ccSums_[i][1]] 
                           << " R=" << sign_*rightConstResults[ccSums_[i][2]]);
    }

  /** Do constant-vector sums */
  for (unsigned int i=0; i<cvSums_.size(); i++)
    {
      Tabs tab1;
      RefCountPtr<EvalVector>& v = rightVectorResults[cvSums_[i][2]];
      SUNDANCE_VERB_MEDIUM(tab1 << "doing c-v sum for " 
                           << vectorResultDeriv(cvSums_[i][0])
                           << ": L=" << leftConstResults[cvSums_[i][1]] 
                           << " R=" << sign_ << "*" 
                           << rightVectorResults[cvSums_[i][2]]->str());
      if (isOne(sign_))
        {
          v->add_S(leftConstResults[cvSums_[i][1]]);
        }
      else
        {
          v->multiply_S_add_S(sign_, leftConstResults[cvSums_[i][1]]);
        }
      vectorResults[cvSums_[i][0]] = v;
    }

  /* Do vector-constant sums */
  for (unsigned int i=0; i<vcSums_.size(); i++)
    {
      Tabs tab1;
      RefCountPtr<EvalVector>& v = leftVectorResults[vcSums_[i][1]] ;
      SUNDANCE_VERB_MEDIUM(tab1 << "doing v-c sum for " 
                           << vectorResultDeriv(vcSums_[i][0])
                           << ": L=" << leftVectorResults[vcSums_[i][1]]->str()
                           << " R=" 
                           << sign_*rightConstResults[vcSums_[i][2]]);
      v->add_S(sign_*rightConstResults[vcSums_[i][2]]);
      vectorResults[vcSums_[i][0]] = v;
    }

  /* Do vector-vector sums */
  for (unsigned int i=0; i<vvSums_.size(); i++)
    {
      Tabs tab1;
      RefCountPtr<EvalVector>& v = leftVectorResults[vvSums_[i][1]];
      SUNDANCE_VERB_MEDIUM(tab1 << "doing v-v sum for " 
                           << vectorResultDeriv(vvSums_[i][0])
                           << ": L=" 
                           << leftVectorResults[vvSums_[i][1]]->str() 
                           << " R=" << sign_ << "*" 
                           << rightVectorResults[vvSums_[i][2]]->str());
      if (isOne(sign_))
        {
          v->add_V(rightVectorResults[vvSums_[i][2]].get());
        }
      else
        {
          v->add_SV(sign_, rightVectorResults[vvSums_[i][2]].get());
        }
      vectorResults[vvSums_[i][0]] = v;
    }
  
  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "sum result " << endl;
      sparsity()->print(cerr, vectorResults,
                        constantResults);
    }
}

