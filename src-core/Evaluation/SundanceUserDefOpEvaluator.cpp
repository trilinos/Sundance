/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUserDefOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUserDefOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;




UserDefOpEvaluator
::UserDefOpEvaluator(const UserDefOp* expr,
                     const EvalContext& context)
  : SubtypeEvaluator<UserDefOp>(expr, context),
    childExpr_(expr->numChildren()),
    childSparsity_(expr->numChildren()),
    childEval_(expr->numChildren()),
    maxOrder_(0),
    d0ResultIndex_(-1),
    d0ArgDerivIndex_(expr->numChildren()),
    d0ArgDerivIsConstant_(expr->numChildren()),
    constantArgPtr_(),
    vectorArgPtr_()
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing user defined op evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  for (int i=0; i<expr->numChildren(); i++)
    {
      childExpr_[i] = expr->evaluatableChild(i);
      childSparsity_[i] = childExpr_[i]->sparsitySuperset(context);
      childEval_[i] = childExpr_[i]->evaluator(context);
      childEval_[i]->addClient();
      SUNDANCE_VERB_MEDIUM(tabs << "argument #" << i << " sparsity " << endl 
                           << *(childSparsity_[i]));
    }

  for (unsigned int c=0; c<childSparsity_.size(); c++)
    {
      /* Find the index of the argument's value (zeroth-order deriv).
       * If this does not exist, we have an error.  
       */
      
      d0ArgDerivIndex_[c] = -1;
      const RefCountPtr<SparsitySuperset>& argSparsity = childSparsity_[c];
      for (int i=0; i<argSparsity->numDerivs(); i++)
        {
          const MultipleDeriv& d = argSparsity->deriv(i);
          
          if (d.order()==0)
            {
              d0ArgDerivIndex_[c] = i;
              if (argSparsity->state(i)==ConstantDeriv)
                {
                  d0ArgDerivIsConstant_[c] = true;
                  constantArgPtr_.append(c);
                }
              else
                {
                  d0ArgDerivIsConstant_[c] = false;
                  vectorArgPtr_.append(c);
                }
            }
        }
      
      TEST_FOR_EXCEPTION(d0ArgDerivIndex_[c] == -1, InternalError,
                         "Inconsistency in sparsity patterns of nonlin "
                         "operator " + expr->toString() + " and its "
                         "argument. "
                         "The zeroth-order derivative of the argument "
                         "was not found in the sparsity pattern "
                         + argSparsity->toString());
    }
         
         
  /*
   * We now build up the tables specifying how to compute the
   * various derivatives.
   */
  
  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity()->deriv(i);
      TEST_FOR_EXCEPTION(d.order() > 0, RuntimeError,
                         "deriv order > 0 not implemented for "
                         "user-defined ops");
      
      /* We need to keep track of the max order, so we know how many
       * derivatives of the operand to compute */
      
      int order = d.order();
      
      if (order > maxOrder_) maxOrder_ = order;
      
      if (order==0)
        {
          /* The zeroth order functional derivative is just the operator
           * F(f,g,h,...), so we need not record any derivatives of the 
           * argument. All we have to do is record the index into which
           * we will write the result of the operation. */
          
          d0ResultIndex_ = i;
          addVectorIndex(i, i); /* result is a vector */
        }
    }

  TEST_FOR_EXCEPTION(d0ResultIndex_==-1, InternalError,
                     "deriv order zero not found in user def op evaluator"
                     "ctor");
                   

}

void UserDefOpEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  TimeMonitor timer(evalTimer());
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tabs << "------- UserDefOpEvaluator::eval() -------");

  /* evaluate the arguments */
  Array<Array<RefCountPtr<EvalVector> > > argVectorResults(childEval_.size());
  Array<Array<double> > argConstantResults(childEval_.size());

  evalChildren(mgr, argConstantResults, argVectorResults);

  Array<RefCountPtr<EvalVector> > funcVectorArgs;
  Array<double> funcConstantArgs;
  
  int numPoints = 1;
  for (int c=0; c<d0ArgDerivIndex_.size(); c++)
    {
      int valPtr = d0ArgDerivIndex_[c];
      if (d0ArgDerivIsConstant_[c])
        {
          funcConstantArgs.append(argConstantResults[c][valPtr]);
        }
      else
        {
          funcVectorArgs.append(argVectorResults[c][valPtr]);
          numPoints = argVectorResults[c][valPtr]->length();
        }
    }

  RefCountPtr<EvalVector> opDerivs = mgr.popVector();
  evalOperator(numPoints, funcConstantArgs, funcVectorArgs, 
               constantArgPtr_, vectorArgPtr_,
               opDerivs);

  /*
   * Allocate the results array
   */
  vectorResults.resize(sparsity()->numDerivs());
  

  /* --- Zeroth derivative term */

  vectorResults[d0ResultIndex_] = opDerivs;
  
}

void UserDefOpEvaluator::resetNumCalls() const
{
  for (unsigned int i=0; i<childEval_.size(); i++)
    {
      childEval_[i]->resetNumCalls();
    }
  Evaluator::resetNumCalls();
}

void UserDefOpEvaluator::evalChildren(const EvalManager& mgr,
                                      Array<Array<double> >& constResults,
                                      Array<Array<RefCountPtr<EvalVector> > >& vecResults) const 
  
{
  for (unsigned int i=0; i<childEval_.size(); i++)
    {
      childEval_[i]->eval(mgr, constResults[i], vecResults[i]);
    }
}



void UserDefOpEvaluator
::evalOperator(int numPoints,
               const Array<double>& constantArg,
               const Array<RefCountPtr<EvalVector> >& vectorArg,
               const Array<int>& constantArgPtr,
               const Array<int>& vectorArgPtr,
               RefCountPtr<EvalVector>& opResults) const
{
  Array<double> evalPt(constantArgPtr.size() + vectorArgPtr.size());

  opResults->resize(numPoints);

  for (int p=0; p<numPoints; p++)
    {
      for (unsigned int c=0; c<constantArgPtr.size(); c++)
        {
          evalPt[constantArgPtr[c]] = constantArg[c];
        }
      for (unsigned int v=0; v<vectorArgPtr.size(); v++)
        {
          evalPt[vectorArgPtr[v]] = vectorArg[v]->start()[p];
        }
      opResults->start()[p] = expr()->op()->eval0(evalPt);
    }
}


