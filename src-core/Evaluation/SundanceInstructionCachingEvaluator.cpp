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

#include "SundanceInstructionCachingEvaluator.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceUnaryMinus.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceNonlinearUnaryOp.hpp"
#include "SundanceUserDefOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

InstructionCachingEvaluatorFactory::InstructionCachingEvaluatorFactory()
  : EvaluatorFactory()
{;}

Evaluator* InstructionCachingEvaluatorFactory
::createEvaluator(const EvaluatableExpr* expr,
                  int derivSetIndex) const
{
  /* do brute-force double dispatch to create 
   * the appropriate evaluator subtype */


  const SumExpr* s = dynamic_cast<const SumExpr*>(expr);
  if (s != 0)
    {
      return new InstructionCachingSumEvaluator(s, derivSetIndex);
    }

  const ProductExpr* p = dynamic_cast<const ProductExpr*>(expr);
  if (p != 0)
    {
      return new InstructionCachingProductEvaluator(p, derivSetIndex);
    }
  
  const DiffOp* d = dynamic_cast<const DiffOp*>(expr);
  if (d != 0)
    {
      return new BruteForceDiffOpEvaluator(d);
    }

  const UnaryMinus* um = dynamic_cast<const UnaryMinus*>(expr);
  if (um != 0)
    {
      return new InstructionCachingUnaryMinusEvaluator(um, derivSetIndex);
    }

  const NonlinearUnaryOp* ue = dynamic_cast<const NonlinearUnaryOp*>(expr);
  if (ue != 0)
    {
      return new InstructionCachingNonlinearUnaryOpEvaluator(ue, derivSetIndex);
    }


  const UserDefOp* ude = dynamic_cast<const UserDefOp*>(expr);
  if (ude != 0)
    {
      return new BruteForceUserDefOpEvaluator(ude);
    }
  /** The expr seems to be one that can be handled by the
   * base class, so forward to the base class' factory method */
  return commonCreate(expr);

  TEST_FOR_EXCEPTION(true, InternalError,
                     "InstructionCachingEvaluatorFactory could not create an "
                     "evaluator for " << expr->toString());

  return new InstructionCachingSumEvaluator(s, derivSetIndex); // return dummy to placate -Wall 
  
}


InstructionCachingSumEvaluator
::InstructionCachingSumEvaluator(const SumExpr* expr,
                                 int derivSetIndex)
: SumEvaluator(expr), 
  leftIndex_(), 
  rightIndex_(), 
  derivSetIndex_(derivSetIndex),
  leftDerivSetIndex_(), 
  rightDerivSetIndex_(),
  isNegative_(expr->sign() < 0), 
  leftExpr_(expr->leftEvaluatable()), 
  rightExpr_(expr->rightEvaluatable()), 
  sparsity_(expr->sparsity(derivSetIndex).get()),
  needsInit_(true)
{
  verbosity() = Evaluator::classVerbosity();
}

void InstructionCachingSumEvaluator::init()
{
  Tabs tabs;
  leftDerivSetIndex_ = expr()->leftDerivSetIndex(derivSetIndex_);
  rightDerivSetIndex_ = expr()->rightDerivSetIndex(derivSetIndex_);
  
  const SparsityPattern* leftSparsity 
    = leftExpr_->sparsity(leftDerivSetIndex_).get();

  const SparsityPattern* rightSparsity 
    = rightExpr_->sparsity(rightDerivSetIndex_).get();

  if (verbosity() > VerbLow) 
    {
      Tabs tabs;
      cerr << tabs << "initializing IC sum evaluator" << endl;
      cerr << tabs << "left deriv set index = " << leftDerivSetIndex_ << endl;
      cerr << tabs << "left sparsity" << endl << *leftSparsity << endl;
      cerr << tabs << "right deriv set index = " 
           << rightDerivSetIndex_ << endl;
      cerr << tabs << "right sparsity" << endl << *rightSparsity << endl;
    }

  for (int i=0; i<sparsity_->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity_->deriv(i);
      TEST_FOR_EXCEPTION(!leftSparsity->containsDeriv(d) 
                         && !rightSparsity->containsDeriv(d),
                         InternalError,
                         "deriv " << d.toString() 
                         << " was not found in either left or right operand "
                         "of expr " << expr()->toString());
      
      int iLeft = -1;
      int iRight = -1;
      if (leftSparsity->containsDeriv(d)) iLeft = leftSparsity->getIndex(d);
      if (rightSparsity->containsDeriv(d)) iRight = rightSparsity->getIndex(d);

      leftIndex_.append(iLeft);
      rightIndex_.append(iRight);
    }


  if (verbosity() > VerbSilent) 
    {
      cerr << tabs << "instruction tables: " << endl;
      cerr << tabs << "left: " << leftIndex_ << endl;
      cerr << tabs << "right: " << rightIndex_ << endl;
    }

}

void InstructionCachingSumEvaluator
::eval(const EvalManager& mgr,
       RefCountPtr<EvalVectorArray>& results) const
{ 
  TimeMonitor timer(sumEvalTimer());
  Tabs tabs;


  SUNDANCE_OUT(this->verbosity() > VerbLow,
               tabs << "------- InstructionCachingSumEvaluator -------");

  if (needsInit_) 
    {
      const_cast<InstructionCachingSumEvaluator*>(this)->init(); 
      needsInit_=false;
    }


  results = mgr.stack().popVectorArray(sparsity_);

  RefCountPtr<EvalVectorArray> leftResults; 
  RefCountPtr<EvalVectorArray> rightResults; 

  SUNDANCE_OUT(this->verbosity() > VerbLow, 
               tabs << "Evaluating left");
  leftExpr_->evaluate(mgr, leftResults);

  SUNDANCE_OUT(this->verbosity() > VerbLow, 
               tabs << "Evaluating right");
  rightExpr_->evaluate(mgr, rightResults);

  for (int i=0; i<results->size(); i++)
    {
      int iLeft = leftIndex_[i];
      int iRight = rightIndex_[i];

      if (iLeft == -1)
        {
          (*results)[i]->copy((*rightResults)[iRight]);
          if (isNegative_ && !(*results)[i]->isZero()) 
            {
              (*results)[i]->unaryMinus();
            }
        }
      else if (iRight == -1)
        {
          (*results)[i]->copy((*leftResults)[iLeft]);
        }
      else
        {
          (*results)[i]->copy((*leftResults)[iLeft]);
          if (!isNegative_) 
            {
              (*results)[i]->addScaled((*rightResults)[iRight], 1.0);
            }
          else
            { 
              (*results)[i]->addScaled((*rightResults)[iRight], -1.0);
            }
        }
    }

  if (verbosity() > VerbLow)
    {
      cerr << endl << tabs << "sum eval results: " << endl;
      {
        Tabs t2;
        cerr << endl << t2 << "left results: " << endl;
        leftResults->print(cerr, 
                           leftExpr_->getDerivSet(leftDerivSetIndex_));
        cerr << endl << t2 << "right results: " << endl;
        rightResults->print(cerr, 
                            rightExpr_->getDerivSet(rightDerivSetIndex_));
      }
      cerr << endl << tabs << "final results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex_));
    }
}




InstructionCachingProductEvaluator
::InstructionCachingProductEvaluator(const ProductExpr* expr,
                                     int derivSetIndex)
  : ProductEvaluator(expr), 
  leftIndex_(), 
  rightIndex_(), 
  derivSetIndex_(derivSetIndex),
  leftDerivSetIndex_(), 
  rightDerivSetIndex_(),
  leftExpr_(expr->leftEvaluatable()), 
  rightExpr_(expr->rightEvaluatable()), 
  sparsity_(expr->sparsity(derivSetIndex).get()),
  needsInit_(true)
{
  verbosity() = Evaluator::classVerbosity();
}

void InstructionCachingProductEvaluator::init()
{
  Tabs tabs;

  leftDerivSetIndex_ = expr()->leftDerivSetIndex(derivSetIndex_);
  rightDerivSetIndex_ = expr()->rightDerivSetIndex(derivSetIndex_);

  const SparsityPattern* leftSparsity 
    = expr()->leftEvaluatable()->sparsity(leftDerivSetIndex_).get();

  const SparsityPattern* rightSparsity 
    = expr()->rightEvaluatable()->sparsity(rightDerivSetIndex_).get();

  if (verbosity() > VerbLow) 
    {
      cerr << tabs << "initializing IC product evaluator" << endl;
      cerr << tabs << "left sparsity" << endl << *leftSparsity << endl;
      cerr << tabs << "right sparsity" << endl << *rightSparsity << endl;
    }

  for (int i=0; i<sparsity_->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity_->deriv(i);
      
      if (d.order()==0)
        {
          TEST_FOR_EXCEPTION(!leftSparsity->containsDeriv(d),
                             InternalError,
                             "InstructionCachingProductEvaluator ctor: "
                             "derivative " << d << " not found in "
                             "left sparsity pattern " << *leftSparsity);
          TEST_FOR_EXCEPTION(!rightSparsity->containsDeriv(d),
                             InternalError,
                             "InstructionCachingProductEvaluator ctor: "
                             "derivative " << d << " not found in "
                             "right sparsity pattern " << *rightSparsity);
          
          int iLeft = leftSparsity->getIndex(d);
          int iRight = rightSparsity->getIndex(d);
          SUNDANCE_OUT(this->verbosity() > VerbLow, 
                       "indices of left and right results vectors: L="
                       << iLeft << " R=" << iRight);
          leftIndex_.append(tuple(iLeft));
          rightIndex_.append(tuple(iRight));
      
      
        }
      else
        {
          Array<MultipleDeriv> leftOps;
          Array<MultipleDeriv> rightOps;

          d.productRulePermutations(leftOps, rightOps);
          
          Tabs t1;

          Array<int> tmpLeft;
          Array<int> tmpRight;
          for (int j=0; j<leftOps.size(); j++)
            {
              Tabs t2;
              
              const MultipleDeriv& dLeft = leftOps[j];
              const MultipleDeriv& dRight = rightOps[j];
              if (!leftSparsity->containsDeriv(dLeft)
                  || !rightSparsity->containsDeriv(dRight)) continue;

              int iLeft = leftSparsity->getIndex(dLeft);
              int iRight = rightSparsity->getIndex(dRight);
              tmpLeft.append(iLeft);
              tmpRight.append(iRight);
            }
          leftIndex_.append(tmpLeft);
          rightIndex_.append(tmpRight);
        }
    }


  if (verbosity() > VerbSilent) 
    {
      cerr << tabs << "instruction tables: " << endl;
      cerr << tabs << "left: " << leftIndex_ << endl;
      cerr << tabs << "right: " << rightIndex_ << endl;
    }
  
}

void InstructionCachingProductEvaluator
::eval(const EvalManager& mgr,
       RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(productEvalTimer());
  Tabs tabs;

  SUNDANCE_OUT(this->verbosity() > VerbLow,
               tabs << "------- InstructionCachingProductEvaluator -------");


  if (needsInit_) 
    {
      const_cast<InstructionCachingProductEvaluator*>(this)->init(); 
      needsInit_=false;
    }

  results = mgr.stack().popVectorArray(sparsity_);

  RefCountPtr<EvalVectorArray> leftResults; 
  RefCountPtr<EvalVectorArray> rightResults; 

  if (verbosity() > VerbLow)
    {
      cerr << endl << tabs << "eval left operand" << endl;
    }
  leftExpr_->evaluate(mgr, leftResults);

  if (verbosity() > VerbLow)
    {
      cerr << endl << tabs << "eval right operand" << endl;
    }
  rightExpr_->evaluate(mgr, rightResults);

  if (verbosity() > VerbMedium)
    {
      cerr << endl << tabs << "product operand results: " << endl;
      {
        Tabs t2;
        cerr << endl << t2 << "left results: " << endl;
        leftResults->print(cerr, 
                           leftExpr_->getDerivSet(leftDerivSetIndex_));
        cerr << endl << t2 << "right results: " << endl;
        rightResults->print(cerr, 
                            rightExpr_->getDerivSet(rightDerivSetIndex_));
      }
    }

  for (int i=0; i<results->size(); i++)
    {
      for (int j=0; j<leftIndex_[i].size(); j++)
        {
          int iLeft = leftIndex_[i][j];
          int iRight = rightIndex_[i][j];
          if (j==0)
            {
              (*results)[i]->copy((*leftResults)[iLeft]);
              (*results)[i]->multiply((*rightResults)[iRight]);
            }
          else
            {
              (*results)[i]->addProduct( (*leftResults)[iLeft],
                                         (*rightResults)[iRight]);
            }
        }
    }
     
  if (verbosity() > VerbMedium)
    {
      cerr << endl << tabs << "product eval results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex_));
    }
}



InstructionCachingUnaryMinusEvaluator
::InstructionCachingUnaryMinusEvaluator(const UnaryMinus* expr,
                                        int derivSetIndex)
  : UnaryMinusEvaluator(expr), derivSetIndex_(derivSetIndex),
    sparsity_(expr->sparsity(derivSetIndex).get())
{
  verbosity() = Evaluator::classVerbosity();
}

void InstructionCachingUnaryMinusEvaluator::eval(const EvalManager& mgr,
                                         RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(unaryMinusEvalTimer());
  Tabs tab;
  SUNDANCE_OUT(this->verbosity() > VerbLow,
               tab << "------- InstructionCachingUnaryMinusEvaluator -------");

  results = mgr.stack().popVectorArray(sparsity_);

  RefCountPtr<EvalVectorArray> argResults; 

  SUNDANCE_OUT(this->verbosity() > VerbMedium,
               tab << "eval operand");
  expr()->evaluatableArg()->evaluate(mgr, argResults);

  TEST_FOR_EXCEPTION(results->size() != argResults->size(),
                     InternalError,
                     "InstructionCachingUnaryMinusEvaluator::eval(): "
                     "output results size="
                     << results->size() 
                     << " is not equal to operand results size="
                     << argResults->size());

  results->copy(argResults);

  for (int i=0; i<results->size(); i++)
    {
      (*results)[i]->unaryMinus();
    }
}


InstructionCachingNonlinearUnaryOpEvaluator
::InstructionCachingNonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr,
                                              int derivSetIndex)
  : NonlinearUnaryOpEvaluator(expr), derivSetIndex_(derivSetIndex),
    sparsity_(expr->sparsity(derivSetIndex).get())
{
  verbosity() = Evaluator::classVerbosity();
  maxOrder_ = 0;
  zeroDerivIndex_ = -1;
  for (int i=0; i<sparsity_->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity_->deriv(i);
      if (d.order() == 0) zeroDerivIndex_ = i;
      TEST_FOR_EXCEPTION(d.order() > 1, RuntimeError,
                         "deriv order > 1 not implemented for unary math ops");
      int order = d.order();
      if (maxOrder_ < order) maxOrder_ = order;
    }
}

void InstructionCachingNonlinearUnaryOpEvaluator::eval(const EvalManager& mgr,
                                         RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(nonlinearUnaryExprEvalTimer());
  Tabs tab;
  SUNDANCE_OUT(this->verbosity() > VerbLow,
               tab << "------- InstructionCachingNonlinearUnaryOpEvaluator -------");

  results = mgr.stack().popVectorArray(sparsity_);

  RefCountPtr<EvalVectorArray> argResults; 

  SUNDANCE_OUT(this->verbosity() > VerbMedium,
               tab << "eval operand");

  expr()->evaluatableArg()->evaluate(mgr, argResults);

  TEST_FOR_EXCEPTION(results->size() != argResults->size(),
                     InternalError,
                     "InstructionCachingNonlinearUnaryOpEvaluator::eval(): "
                     "output results size="
                     << results->size() 
                     << " is not equal to operand results size="
                     << argResults->size());



  /* Do the evaluation, checking for errors */  
  Array<RefCountPtr<EvalVector> > funcDerivs(maxOrder_+1);
  errno = 0;
  (*argResults)[zeroDerivIndex_]->applyUnaryFunction(expr()->op(),
                                                     funcDerivs);
  TEST_FOR_EXCEPTION(errno==EDOM, RuntimeError,
                     "Domain error in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno==ERANGE, RuntimeError,
                     "Range error in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno==EOVERFLOW, RuntimeError,
                     "Overflow in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno != 0, RuntimeError,
                     "Error " << errno << " in expression " << expr()->toString());
                     

  for (int i=0; i<results->size(); i++)
    {
      const MultipleDeriv& d = sparsity_->deriv(i);
      if (d.order()==0)
        {
          (*results)[i]->copy(funcDerivs[0]);
        }
      else
        {
          (*results)[i]->copy(funcDerivs[1]);
          (*results)[i]->multiply((*argResults)[i]);
        }
    }
 if (verbosity() > 1)
    {
      Tabs tabs;
      cerr << endl << tabs << "nonlinear op eval results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex_));
    }

}
