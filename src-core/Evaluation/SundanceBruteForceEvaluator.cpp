/* @HEADER@ */
/* @HEADER@ */

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

static Time& createEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("evaluator creation"); 
  return *rtn;
}

Time& sumEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("sum eval"); 
  return *rtn;
}

Time& productEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("product eval"); 
  return *rtn;
}

Time& diffOpEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("diff op eval"); 
  return *rtn;
}

Time& unaryMinusEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("unary minus eval"); 
  return *rtn;
}

Time& diffOpCalculationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("diff op calculations"); 
  return *rtn;
}

Time& nonlinearUnaryExprEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("nonlinear unary expr eval"); 
  return *rtn;
}

Time& userDefOpEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("user def expr eval"); 
  return *rtn;
}

BruteForceEvaluatorFactory::BruteForceEvaluatorFactory()
  : EvaluatorFactory()
{;}

Evaluator* BruteForceEvaluatorFactory
::createEvaluator(const EvaluatableExpr* expr,
                  int derivSetIndex) const
{
  TimeMonitor t(createEvalTimer());

  /* do brute-force double dispatch to create 
   * the appropriate evaluator subtype */

  const SumExpr* s = dynamic_cast<const SumExpr*>(expr);
  if (s != 0)
    {
      return new BruteForceSumEvaluator(s);
    }

  const ProductExpr* p = dynamic_cast<const ProductExpr*>(expr);
  if (p != 0)
    {
      return new BruteForceProductEvaluator(p);
    }
  
  const DiffOp* d = dynamic_cast<const DiffOp*>(expr);
  if (d != 0)
    {
      return new BruteForceDiffOpEvaluator(d);
    }

  const UnaryMinus* um = dynamic_cast<const UnaryMinus*>(expr);
  if (um != 0)
    {
      return new BruteForceUnaryMinusEvaluator(um);
    }

  const NonlinearUnaryOp* ue = dynamic_cast<const NonlinearUnaryOp*>(expr);
  if (ue != 0)
    {
      return new BruteForceNonlinearUnaryOpEvaluator(ue);
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
                     "BruteForceEvaluatorFactory could not create an "
                     "evaluator for " << expr->toString());

  return new BruteForceSumEvaluator(s); // return dummy to placate -Wall 
  
}

void BruteForceSumEvaluator::eval(const EvalManager& mgr,
                                  RefCountPtr<EvalVectorArray>& results) const
{ 
  TimeMonitor timer(sumEvalTimer());
  Tabs tabs;

  if (verbosity() > 1) 
    {
      cerr << tabs << "--- BruteForceSumEvaluator --- " << endl;
      cerr << tabs << "left: " << expr()->leftEvaluatable()->toString() 
           << endl;
      cerr << tabs << "right: " << expr()->rightEvaluatable()->toString() 
           << endl;
    }


  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();


  if (verbosity() > 1)
    {
      cerr << tabs << "deriv set index = " << derivSetIndex << endl;
      cerr << tabs << "sparsity = " << endl << *sparsity << endl;
    }
  int leftDerivSetIndex = expr()->leftDerivSetIndex(derivSetIndex);
  int rightDerivSetIndex = expr()->rightDerivSetIndex(derivSetIndex);

  const SparsityPattern* leftSparsity 
    = expr()->leftEvaluatable()->sparsity(leftDerivSetIndex).get();

  const SparsityPattern* rightSparsity 
    = expr()->rightEvaluatable()->sparsity(rightDerivSetIndex).get();

  if (verbosity() > 1) 
    {
      cerr << tabs << "left deriv set index = " << leftDerivSetIndex << endl;
      cerr << tabs << "left sparsity" << endl << *leftSparsity << endl;
      cerr << tabs << "right deriv set index = " << rightDerivSetIndex << endl;
      cerr << tabs << "right sparsity" << endl << *rightSparsity << endl;
    }

  results = mgr.stack().popVectorArray(sparsity);

  RefCountPtr<EvalVectorArray> leftResults; 
  RefCountPtr<EvalVectorArray> rightResults; 

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tabs << "Evaluating left");
  expr()->leftEvaluatable()->evaluate(mgr, leftResults);

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tabs << "Evaluating right");
  expr()->rightEvaluatable()->evaluate(mgr, rightResults);

  for (int i=0; i<results->size(); i++)
    {
      const MultipleDeriv& d = sparsity->deriv(i);
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

      if (iLeft == -1)
        {
          (*results)[i]->copy((*rightResults)[iRight]);
          if (expr()->sign()<0 && !(*results)[i]->isZero()) (*results)[i]->unaryMinus();
        }
      else if (iRight == -1)
        {
          (*results)[i]->copy((*leftResults)[iLeft]);
        }
      else
        {
          (*results)[i]->copy((*leftResults)[iLeft]);
          if (verbosity() > 1) cerr << endl 
                                    << tabs << "sign=" 
                                    << expr()->sign() << endl;
          if (expr()->sign() > 0) 
            {

              (*results)[i]->addScaled((*rightResults)[iRight], 1.0);
            }
          else
            { 
              (*results)[i]->addScaled((*rightResults)[iRight], -1.0);
            }
        }
    }

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "sum eval results: " << endl;
      {
        Tabs t2;
        cerr << endl << t2 << "left results: " << endl;
        leftResults->print(cerr, 
                           expr()->leftEvaluatable()->getDerivSet(leftDerivSetIndex));
        cerr << endl << t2 << "right results: " << endl;
        rightResults->print(cerr, 
                            expr()->rightEvaluatable()->getDerivSet(rightDerivSetIndex));
      }
      cerr << endl << tabs << "final results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex));
    }
}

void BruteForceProductEvaluator::eval(const EvalManager& mgr,
                                      RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(productEvalTimer());
  Tabs tabs;

  if (verbosity() > 1) 
    {
      cerr << tabs << "--- BruteForceProductEvaluator --- " << endl;
      cerr << tabs << "left: " << expr()->leftEvaluatable()->toString() 
           << endl;
      cerr << tabs << "right: " << expr()->rightEvaluatable()->toString() 
           << endl;
    }

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();

  int leftDerivSetIndex = expr()->leftDerivSetIndex(derivSetIndex);
  int rightDerivSetIndex = expr()->rightDerivSetIndex(derivSetIndex);

  const SparsityPattern* leftSparsity 
    = expr()->leftEvaluatable()->sparsity(leftDerivSetIndex).get();

  const SparsityPattern* rightSparsity 
    = expr()->rightEvaluatable()->sparsity(rightDerivSetIndex).get();

  if (verbosity() > 1) 
    {
      cerr << tabs << "left sparsity" << endl << *leftSparsity << endl;
      cerr << tabs << "right sparsity" << endl << *rightSparsity << endl;
    }

  results = mgr.stack().popVectorArray(sparsity);

  RefCountPtr<EvalVectorArray> leftResults; 
  RefCountPtr<EvalVectorArray> rightResults; 

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "eval left operand" << endl;
    }
  expr()->leftEvaluatable()->evaluate(mgr, leftResults);

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "eval right operand" << endl;
    }
  expr()->rightEvaluatable()->evaluate(mgr, rightResults);

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "product operand results: " << endl;
      {
        Tabs t2;
        cerr << endl << t2 << "left results: " << endl;
        leftResults->print(cerr, 
                           expr()->leftEvaluatable()->getDerivSet(leftDerivSetIndex));
        cerr << endl << t2 << "right results: " << endl;
        rightResults->print(cerr, 
                            expr()->rightEvaluatable()->getDerivSet(rightDerivSetIndex));
      }
    }

  for (int i=0; i<results->size(); i++)
    {
      const MultipleDeriv& d = sparsity->deriv(i);

      if (d.order()==0)
        {
          TEST_FOR_EXCEPTION(!leftSparsity->containsDeriv(d),
                             InternalError,
                             "BruteForceProductEvaluator::eval(): "
                             "derivative " << d << " not found in "
                             "left sparsity pattern " << *leftSparsity);
          TEST_FOR_EXCEPTION(!rightSparsity->containsDeriv(d),
                             InternalError,
                             "BruteForceProductEvaluator::eval(): "
                             "derivative " << d << " not found in "
                             "right sparsity pattern " << *rightSparsity);

          int iLeft = leftSparsity->getIndex(d);
          int iRight = rightSparsity->getIndex(d);
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       "indices of left and right results vectors: L="
                       << iLeft << " R=" << iRight);
          
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tabs << "d=0 left=" 
                       << (*leftResults)[iLeft]->getStringValue()
                       << " right=" 
                       << (*rightResults)[iRight]->getStringValue());
                                                    
          (*results)[i]->copy((*leftResults)[iLeft]);
          (*results)[i]->multiply((*rightResults)[iRight]);
        }
      else
        {
          Array<MultipleDeriv> leftOps;
          Array<MultipleDeriv> rightOps;

          d.productRulePermutations(leftOps, rightOps);
          
          (*results)[i]->setToZero();

          Tabs t1;

          if (verbosity() > 1)
            {
              cerr << t1 << "doing product rule sum" << endl;
            }
          for (int j=0; j<leftOps.size(); j++)
            {
              Tabs t2;
              
              const MultipleDeriv& dLeft = leftOps[j];
              const MultipleDeriv& dRight = rightOps[j];
              if (!leftSparsity->containsDeriv(dLeft)
                  || !rightSparsity->containsDeriv(dRight)) continue;

              int iLeft = leftSparsity->getIndex(dLeft);
              int iRight = rightSparsity->getIndex(dRight);
              
              if (verbosity() > 1)
                {
                  cerr << t2 << "left deriv = " << dLeft
                       << "=" << (*leftResults)[iLeft]->getStringValue()
                       << " right deriv = " << dRight
                       << "=" << (*rightResults)[iRight]->getStringValue()
                       << endl;
                }
              (*results)[i]->addProduct( (*leftResults)[iLeft],
                                         (*rightResults)[iRight]);
              if (verbosity() > 1)
                {
                  cerr << t2 << "product rule partial sum = "
                       << (*results)[i]->getStringValue() << endl;
                }
            }
          if (verbosity() > 1)
            {
              cerr << t1 << "done product rule sum: value = " 
                   << (*results)[i]->getStringValue() << endl;
            }
        }
    }

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "product eval results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex));
    }
}

void BruteForceDiffOpEvaluator::eval(const EvalManager& mgr,
                                     RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(diffOpEvalTimer());
  static RefCountPtr<Time> t0 
    = TimeMonitor::getNewTimer("diff op set to zero"); 
  Tabs tabs;

  if (verbosity() > 1) 
    {
      cerr << tabs << "--- BruteForceDiffOpEvaluator --- " << endl;
      cerr << tabs << "op: " << expr()->mi().toString()
           << endl;
      cerr << tabs << "arg: " << expr()->evaluatableArg()->toString() 
           << endl;
    }

  /* look up the sparsity pattern for the result */
  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();

  if (verbosity() > 1) 
    {
      cerr << tabs << "diff op result sparsity" << endl << *sparsity << endl;
    }


  /* look up the sparsity pattern for the argument */
  int argDerivSetIndex = expr()->argDerivSetIndex(derivSetIndex);
  const SparsityPattern* argSparsity 
    = expr()->evaluatableArg()->sparsity(argDerivSetIndex).get();

  if (verbosity() > 1) 
    {
      cerr << tabs << "diff op arg sparsity" << endl << *argSparsity << endl;
    }
  results = mgr.stack().popVectorArray(sparsity);

  RefCountPtr<EvalVectorArray> argResults; 


  /* evaluate the argument */
  if (verbosity() > 1) 
    {
      cerr << tabs << "evaluating arg " << endl;
    }

  {
    static RefCountPtr<Time> t 
      = TimeMonitor::getNewTimer("diff op arg eval"); 
    TimeMonitor timer2(*t);
    expr()->evaluatableArg()->evaluate(mgr, argResults);
  }

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "arg results: " << endl;
      argResults->print(cerr, expr()->evaluatableArg()->getDerivSet(argDerivSetIndex));
    }


  SundanceUtils::Set<Deriv> argDeps;
  expr()->evaluatableArg()->getRoughDependencies(argDeps);

  /* first assemble a list of all functions that we have to evaluate
   * in order to do this diff op */
  SundanceUtils::Set<Deriv> requiredFuncs;
  for (int i=0; i<results->size(); i++)
    {
      if (sparsity->isZero(i))
        {
          continue;
        }
      const MultipleDeriv& d = sparsity->deriv(i);

      if (!expr()->requiresFunctionsToEval(d)) continue;

      const SundanceUtils::Set<Deriv>& funcsToDoD = expr()->requiredFunctions(d);
      SundanceUtils::Set<Deriv>::const_iterator iter;
      for (iter=funcsToDoD.begin(); iter != funcsToDoD.end(); iter++)
        {
          requiredFuncs.put(*iter);
        }
    }

  Map<Deriv, int> derivToFuncResultIndexMap;
  Array<RefCountPtr<EvalVector> > funcResults;
  SundanceUtils::Set<Deriv>::const_iterator funcIter;

  int funcCounter = 0;

  for (funcIter=requiredFuncs.begin(); funcIter != requiredFuncs.end();
       funcIter++)
    {
      const Deriv& f = *funcIter;
      TEST_FOR_EXCEPTION(!f.isFunctionalDeriv(), InternalError,
                         "BruteForceDiffOpEvaluator::eval() detected coord "
                         "deriv where a functional deriv was expected");
      const FuncElementBase* func = f.funcDeriv()->func();
      const MultiIndex& beta = f.funcDeriv()->multiIndex();
      const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(func);
      RefCountPtr<EvalVector> uResult;
      if (u != 0)
        {
          const DiscreteFuncElement* u0 
            = dynamic_cast<const DiscreteFuncElement*>(u->evalPt());
          if (u0 != 0)
            {
              static RefCountPtr<Time> t 
                = TimeMonitor::getNewTimer("diff op discrete func eval"); 
              TimeMonitor timer2(*t);
              uResult = mgr.stack().popFullVector();
              mgr.evalDiscreteFuncElement(u0, beta, uResult);
            }
          else
            {
              TimeMonitor timer2(*t0);
              uResult = mgr.stack().popTrivialVector();
              uResult->setToZero();
            }
          funcResults.append(uResult);
          derivToFuncResultIndexMap.put(f, funcCounter);
          funcCounter++;
        }
      
    }
         
  if (verbosity() > 1)
    {
      cerr << tabs << "applying functional/spatial chain rule "  << endl;
    }

  for (int i=0; i<results->size(); i++)
    {
      Tabs tab0;
      if (sparsity->isZero(i))
        {
          TimeMonitor timer2(*t0);
          (*results)[i]->setToZero();
        }
      else
        {
          const MultipleDeriv& d = sparsity->deriv(i);
          if (verbosity() > 1)
            {
              cerr << tab0 << "computing deriv " << d << endl;
            }
 
          /* first get the deriv wrt {d, me} */
          MultipleDeriv dAndMe = d;
          dAndMe.put(expr()->myCoordDeriv());

          if (argSparsity->containsDeriv(dAndMe))
            {
              Tabs tab1;
              if (verbosity() > 1)
                {
                  cerr << tab1 << "found nonzero " << dAndMe << endl;
                }
              int iArg = argSparsity->getIndex(dAndMe);

              static RefCountPtr<Time> t 
                = TimeMonitor::getNewTimer("diff op arg result copy"); 
              TimeMonitor timer2(*t);
              (*results)[i]->copy((*argResults)[iArg]);
            }
          else
            {
              Tabs tab1;
              if (verbosity() > 1)
                {
                  cerr << tab1 << "deriv " << dAndMe << " is zero" << endl;
                }
              TimeMonitor timer2(*t0);
              (*results)[i]->setToZero();
            }

          /* now add in the terms from the chain rule */
          Array<MultipleDeriv> leftOps;
          Array<MultipleDeriv> rightOps;

          d.productRulePermutations(leftOps, rightOps);
          SundanceUtils::Set<Deriv>::const_iterator iter;

          for (iter=argDeps.begin(); iter != argDeps.end(); iter++)
            {
              Tabs tab1;

              const Deriv& dep = *iter;
                  
              if (verbosity() > 1)
                {
                  cerr << tab1 << "processing arg deriv " 
                       << dep.toString() << endl;
                }

              /* Paranoiacally check that this dependency is 
               * a functional deriv */
              TEST_FOR_EXCEPTION(!dep.isFunctionalDeriv(), InternalError,
                                 "BruteForceDiffOpEvaluator::eval() "
                                 "detected "
                                 "a non-functional derivative in its "
                                 "dependency set");
                             
              const FunctionalDeriv* fDep = dep.funcDeriv();
              const MultiIndex& beta = fDep->multiIndex();
              int depFuncID = fDep->funcID();

              for (int j=0; j<leftOps.size(); j++)
                {
                  if (rightOps[j].order() > 1) continue;

                  if (rightOps[j].order()==1)
                    {
                      Deriv r = *(rightOps[j].begin());
                      if (!r.isFunctionalDeriv()) continue;
                      const FunctionalDeriv* fr = r.funcDeriv();
                      const MultiIndex& gamma = fr->multiIndex();
                      int funcID = fr->funcID();
                      if (funcID != depFuncID || !(expr()->mi()+beta == gamma)) continue;
                      /* HERE IS WHERE WE NEED D_gamma u(fr) */
                    }
                  
                  MultipleDeriv L = leftOps[j];
                  L.put(dep);

                  if (!argSparsity->containsDeriv(L)) 
                    {
                      if (verbosity() > 1)
                        {
                          Tabs tab2;
                          cerr << tab2 << "skipping nonexistent deriv " 
                               << L << endl;
                        }
                      continue;
                    }

                  int iLeft = argSparsity->getIndex(L);
                  if (rightOps[j].order()==0)
                    {
                      Tabs tab2;
                      Deriv fieldDeriv 
                        = fDep->derivWrtMultiIndex(expr()->mi());
                      if (!derivToFuncResultIndexMap.containsKey(fieldDeriv))
                        continue;
                      int fieldIndex 
                        = derivToFuncResultIndexMap.get(fieldDeriv);
                          
                      if (verbosity() > 1)
                        {
                          cerr << tab2 << "adding product (" 
                               << (*argResults)[iLeft]->getStringValue()
                               << ",  " 
                               << funcResults[fieldIndex]->getStringValue() << ") to " << (*results)[i]->getStringValue() << endl;
                        }
                      static RefCountPtr<Time> t 
                        = TimeMonitor::getNewTimer("diff op add product"); 
                      TimeMonitor timer2(*t);
                      (*results)[i]->addProduct((*argResults)[iLeft],
                                                funcResults[fieldIndex]);
                    }
                  else
                    {
                      Tabs tab2;
                      if (verbosity() > 1)
                        {
                          cerr << tab2 << "adding term " 
                               << (*argResults)[iLeft]->getStringValue()
                               << endl;
                        }
                      static RefCountPtr<Time> t 
                        = TimeMonitor::getNewTimer("diff op add scaled"); 
                      TimeMonitor timer2(*t);
                      (*results)[i]->addScaled((*argResults)[iLeft],
                                               1.0);
                    }
                                                
                }
            }
        }
    } 

  if (verbosity() > 1)
    {
      cerr << endl << tabs << "diff op eval results: " << endl;
      results->print(cerr, expr()->getDerivSet(derivSetIndex));
    }
}

void BruteForceUnaryMinusEvaluator::eval(const EvalManager& mgr,
                                         RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(unaryMinusEvalTimer());
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tab << "------- BruteForceUnaryMinusEvaluator -------");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();

  results = mgr.stack().popVectorArray(sparsity);

  RefCountPtr<EvalVectorArray> argResults; 

  SUNDANCE_OUT(verbosity() > VerbMedium,
               tab << "eval operand");
  expr()->evaluatableArg()->evaluate(mgr, argResults);

  TEST_FOR_EXCEPTION(results->size() != argResults->size(),
                     InternalError,
                     "BruteForceUnaryMinusEvaluator::eval(): "
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


void BruteForceNonlinearUnaryOpEvaluator::eval(const EvalManager& mgr,
                                         RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(nonlinearUnaryExprEvalTimer());
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tab << "------- BruteForceNonlinearUnaryOpEvaluator -------");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();

  results = mgr.stack().popVectorArray(sparsity);

  RefCountPtr<EvalVectorArray> argResults; 

  SUNDANCE_OUT(verbosity() > VerbMedium,
               tab << "eval operand");

  expr()->evaluatableArg()->evaluate(mgr, argResults);

  TEST_FOR_EXCEPTION(results->size() != argResults->size(),
                     InternalError,
                     "BruteForceNonlinearUnaryOpEvaluator::eval(): "
                     "output results size="
                     << results->size() 
                     << " is not equal to operand results size="
                     << argResults->size());



  int maxOrder = 0;
  int zeroDerivIndex = -1;
  for (int i=0; i<results->size(); i++)
    {
      const MultipleDeriv& d = sparsity->deriv(i);
      if (d.order() == 0) zeroDerivIndex = i;
      TEST_FOR_EXCEPTION(d.order() > 1, RuntimeError,
                         "deriv order > 1 not implemented for unary math ops");
      int order = d.order();
      if (maxOrder < order) maxOrder = order;
    }
  
  TEST_FOR_EXCEPTION(zeroDerivIndex < 0, RuntimeError,
                     "no zero-order deriv of argument in unary math op eval");
  

  /* Do the evaluation, checking for errors */  
  Array<RefCountPtr<EvalVector> > funcDerivs(maxOrder+1);
  errno = 0;
  (*argResults)[zeroDerivIndex]->applyUnaryFunction(expr()->op(),
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
      const MultipleDeriv& d = sparsity->deriv(i);
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
      results->print(cerr, expr()->getDerivSet(derivSetIndex));
    }

}



void BruteForceUserDefOpEvaluator::eval(const EvalManager& mgr,
                                        RefCountPtr<EvalVectorArray>& results) const
{
  TimeMonitor timer(userDefOpEvalTimer());
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tab << "------- BruteForceUserDefOpEvaluator -------");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const SparsityPattern* sparsity = expr()->sparsity(derivSetIndex).get();

  results = mgr.stack().popVectorArray(sparsity);

  Array<RefCountPtr<EvalVectorArray> > argResults(expr()->numChildren()); 

  SUNDANCE_OUT(verbosity() > VerbMedium,
               tab << "eval operands");

  for (int i=0; i<expr()->numChildren(); i++)
    {
      expr()->evaluatableChild(i)->evaluate(mgr, argResults[i]);
    }


  int maxOrder = 0;
  int zeroDerivIndex = -1;
  for (int i=0; i<results->size(); i++)
    {
      const MultipleDeriv& d = sparsity->deriv(i);
      if (d.order() == 0) zeroDerivIndex = i;
      TEST_FOR_EXCEPTION(d.order() > 0, RuntimeError,
                         "deriv order > 0 not implemented for user def ops");
      int order = d.order();
      if (maxOrder < order) maxOrder = order;
    }
  
  TEST_FOR_EXCEPTION(zeroDerivIndex < 0, RuntimeError,
                     "no zero-order deriv of argument in unary math op eval");
  
  
  int numPts = 0;
  for (int i=0; i<argResults.size(); i++)
    {
      for (int j=0; j<argResults[i]->size(); j++)
        {
          if (!(*(argResults[i]))[j]->isConstant()) 
            numPts = (*(argResults[i]))[j]->length();
        }
    }

  Array<Array<double> > evalPts(numPts);
  for (int p=0; p<numPts; p++)
    {
      evalPts[p].resize(argResults.size());
      for (int f=0; f<argResults.size(); f++)
        {
          const EvalVectorArray& ar = *(argResults[f]);
          if (ar[0]->isConstant())
            {
              evalPts[p][f] = ar[0]->getConstantValue();
            }
          else
            {
              evalPts[p][f] = (ar[0]->start())[p];
            }
        }
    }
  
  /* Do the evaluation, checking for errors */  
  ((*results)[zeroDerivIndex])->resize(numPts);
  for (int p=0; p<numPts; p++)
    {
      (*results)[zeroDerivIndex]->start()[p] = expr()->op()->eval(evalPts[p]);
    }
  

}
