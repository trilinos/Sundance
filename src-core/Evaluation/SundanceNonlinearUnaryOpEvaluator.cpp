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

#include "SundanceNonlinearUnaryOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceNonlinearUnaryOp.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;




NonlinearUnaryOpEvaluator
::NonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr,
                            const EvalContext& context)
  : UnaryEvaluator<NonlinearUnaryOp>(expr, context),
    maxOrder_(0),
    d0ResultIndex_(-1),
    d0ArgDerivIndex_(-1),
    d0ArgDerivIsConstant_(false),
    d1ResultIndex_(),
    d1ArgDerivIndex_(),
    d1ArgDerivIsConstant_(),
    d2ResultIndex_(),
    d2ArgDerivIndex_(),
    d2ArgDerivIsConstant_()
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing nonlinear unary op evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  SUNDANCE_VERB_MEDIUM(tabs << "argument sparsity subset" << endl 
                       << *(argSparsitySubset()));

  /* Find the index of the argument's value (zeroth-order deriv).
   * If this does not exist, we have an error.  
   */
  
  d0ArgDerivIndex_ = -1;
  for (int i=0; i<argSparsitySubset()->numDerivs(); i++)
    {
      const MultipleDeriv& d = argSparsitySubset()->deriv(i);
      
      if (d.order()==0)
        {
          /* we'll need the index of the argument's value in the superset
           * of results */
          int argIndexInSuperset = argSparsitySuperset()->getIndex(d);
          SUNDANCE_VERB_MEDIUM(tabs << "arg value found at index "
                               << argIndexInSuperset 
                               << "in arg result superset");
          TEST_FOR_EXCEPTION(argIndexInSuperset==-1, InternalError,
                             "derivative " << d
                             << " found in arg subset but not in "
                             "arg superset");
          d0ArgDerivIndex_ = argIndexInSuperset;
          if (argSparsitySubset()->state(i)==ConstantDeriv)
            {
              d0ArgDerivIsConstant_ = true;
              SUNDANCE_VERB_MEDIUM(tabs << "arg value is constant");
            }
          else
            {
              d0ArgDerivIsConstant_ = false;
              SUNDANCE_VERB_MEDIUM(tabs << "arg value is non-constant");
            }
        }
    }

  TEST_FOR_EXCEPTION(d0ArgDerivIndex_ == -1, InternalError,
                     "Inconsistency in sparsity patterns of nonlin "
                     "operator " + expr->toString() + " and its "
                     "argument. "
                     "The zeroth-order derivative of the argument "
                     "was not found in the sparsity pattern "
                     + argSparsitySubset()->toString());


  /*
   * We now build up the tables specifying how to compute the
   * various derivatives.
   */
  
  int constCounter = 0;
  int vecCounter = 0;
  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity()->deriv(i);
      TEST_FOR_EXCEPTION(d.order() > 2, RuntimeError,
                         "deriv order > 2 not implemented for unary math ops");

      /* We need to keep track of the max order, so we know how many
       * derivatives of the operand to compute */

      int order = d.order();

      if (order > maxOrder_) maxOrder_ = order;

      if (order==0)
        {
          /* The zeroth order functional derivative is just the operator
           * F(g), so we need not record any derivatives of the 
           * argument. All we have to do is record the index into which
           * we will write the result of the operation. */
          if (d0ArgDerivIsConstant_)
            {
              d0ResultIndex_ = constCounter;;
              addConstantIndex(i, constCounter++); /* result is a constant */
            }
          else
            {
              d0ResultIndex_ = vecCounter;
              addVectorIndex(i, vecCounter++); /* result is a vector */
            }
          SUNDANCE_VERB_MEDIUM(tabs << "operator value goes into result index"
                               << d0ResultIndex_);
        }
      else if (order==1)
        {
          /* The first order functional derivative wrt u is
           * F_u = F'(g) g_u. We must record the index where we will 
           * write this result, the index into the arg's results
           * from which we obtain the arg's derivative g_u, and whether
           * the arg's deriv g_u is a constant or a vector. 
           *
           * The arg's derivative g_u is not used by any other expression
           * of zeroth or first order (it is used in second derivative
           * evaluation) so provided that we evaluate from
           * higher to lower differentiation, we can safely overwrite
           * that vector, performing the calculation F'(g) g_u in place.
           * This will save us a vector copy. 
           */

          /* Record the index at which we will record this derivative */
          Tabs tab1;
          d1ResultIndex_.append(i);
          addVectorIndex(i, vecCounter++); /* result is a vector */
          
          TEST_FOR_EXCEPTION(!argSparsitySubset()->containsDeriv(d),
                             InternalError,
                             "Inconsistency in sparsity patterns of nonlin "
                             "operator " + expr->toString() + " and its "
                             "argument." 
                             "Operator's derivative " + d.toString() 
                             + " was not found in arg's sparsity pattern ");


          /* Record index of g_u in arg results, and also whether it
           * is a constant or a vector */
          int index = argSparsitySuperset()->getIndex(d);
          SUNDANCE_VERB_MEDIUM(tab1 << "deriv " << d << " found at index "
                               << index << " in arg results");

          if (argSparsitySuperset()->state(index)==ConstantDeriv)
            {
              Tabs tab2;
              SUNDANCE_VERB_MEDIUM(tab2 << "arg deriv is constant");
              TEST_FOR_EXCEPTION(!argEval()->constantIndexMap().containsKey(index),
                                 InternalError,
                                 "index " << index 
                                 << " not found in arg constant index map");
              d1ArgDerivIndex_.append(argEval()->constantIndexMap().get(index));
              d1ArgDerivIsConstant_.append(true);
            }
          else
            {
              Tabs tab2;
              SUNDANCE_VERB_MEDIUM(tab2 << "arg deriv is non-constant");
              TEST_FOR_EXCEPTION(!argEval()->vectorIndexMap().containsKey(index),
                                 InternalError,
                                 "index " << index 
                                 << " not found in arg vector index map");
              d1ArgDerivIndex_.append(argEval()->vectorIndexMap().get(index));
              d1ArgDerivIsConstant_.append(false);
            }
          
        }
      else if (order==2)
        {

          /* A second derivative wrt u and v will look like 
           * F_uv = F"(g) g_u g_v + F' g_uv. 
           * We therefore need to know the indices from which to obtain
           * g_u, g_v, and g_uv, as well as whether each is constant or
           * a full vector. Furthermore, g_uv may be structurally zero,
           * in which case we record an index of -1. 
           * We record the indices of {g_uv, g_u, g_v} (in that order)
           * in the array d2ArgDerivIndex, and a constancy marker
           * for each in d2argDerivIsConstant. 
           *
           * The arg's second derivative g_uv is not used by any other 
           * derivative calculation. We can therefore safely perform the
           * operation F' g_uv + F" g_u g_v in place, overwriting the
           * value of g_uv. Note that because second-order derivs require
           * g_u and g_v, these should always be performed before 
           * first-order calculations so that g_u and g_v can be
           * overwritten as well.
           */
          d2ResultIndex_.append(i);
          addVectorIndex(i, vecCounter++); /* result is a vector */
          

          /* 
           * Find the index of the argument's second deriv, and determine
           * whether it is constant.
           */
          Array<int> d2ArgDerivIndex(3);
          Array<int> d2ArgDerivIsConstant(3);

          if (argSparsitySuperset()->containsDeriv(d))
            {
              int index = argSparsitySuperset()->getIndex(d);
              if (argSparsitySuperset()->state(index)==ConstantDeriv)
                {
                  d2ArgDerivIndex[0] 
                    = argEval()->constantIndexMap().get(index);
                  d2ArgDerivIsConstant[0] = true;
                }
              else
                {
                  d2ArgDerivIndex[0] 
                    = argEval()->vectorIndexMap().get(index);
                  d2ArgDerivIsConstant[0] = false;
                }
            }
          else
            {
              d2ArgDerivIndex[0] = -1;
            }

          /* 
           * Find the indices of the two first-order factors g_u and g_v,
           * and whether thery are constant. 
           */
          int j=0;
          for (MultipleDeriv::const_iterator 
                 iter=d.begin(); iter != d.end(); iter++, j++)
            {
              MultipleDeriv sd;
              sd.put(*iter);
              TEST_FOR_EXCEPTION(!argSparsitySubset()->containsDeriv(sd),
                                 InternalError,
                                 "First deriv " << sd << " not found in "
                                 "sparsity pattern for an argument that "
                                 "contains a second derivative");
              int index = argSparsitySuperset()->getIndex(sd);
              DerivState state = argSparsitySuperset()->state(index);
              if (state == ConstantDeriv)
                {
                  d2ArgDerivIndex[1+j] 
                    = argEval()->constantIndexMap().get(index);
                  d2ArgDerivIsConstant[1+j] = true;
                }
              else
                {
                  d2ArgDerivIndex[1+j] 
                    = argEval()->vectorIndexMap().get(index);
                  d2ArgDerivIsConstant[1+j] = false;
                }
            }
          d2ArgDerivIndex_.append(d2ArgDerivIndex);
          d2ArgDerivIsConstant_.append(d2ArgDerivIsConstant);
        }
    }

}

void NonlinearUnaryOpEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const
{
  TimeMonitor timer(evalTimer());
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow,
               tabs << "NonlinearUnaryOpEvaluator::eval() expr="
               << expr()->toString());

  /* evaluate the argument */
  Array<RefCountPtr<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  evalOperand(mgr, argConstantResults, argVectorResults);




  /* Do the application of the nonlinear operator, checking for errors */  


  Array<RefCountPtr<EvalVector> > opDerivs(maxOrder_+1);
  errno = 0;
  
  /* If the argument is a constant, copy it into a vector. */
  RefCountPtr<EvalVector> argValue;
  if (d0ArgDerivIsConstant_)
    {
      argValue = mgr.popVector();
      argValue->resize(1);
      argValue->setToConstant(argConstantResults[d0ArgDerivIndex_]);
    }
  else
    {
      argValue = argVectorResults[d0ArgDerivIndex_];
    }

  SUNDANCE_VERB_HIGH(tabs << "applying operator");

  argValue->applyUnaryOperator(expr()->op(), opDerivs);

  TEST_FOR_EXCEPTION(errno==EDOM, RuntimeError,
                     "Domain error in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno==ERANGE, RuntimeError,
                     "Range error in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno==EOVERFLOW, RuntimeError,
                     "Overflow in expression " << expr()->toString());
  TEST_FOR_EXCEPTION(errno != 0, RuntimeError,
                     "Error " << errno << " in expression " << expr()->toString());

  if (d0ArgDerivIsConstant_)
    {
      constantResults.resize(1);
      constantResults[0] = opDerivs[0]->start()[0];
    }
  else
    {
  
      /*
       * Allocate the results array
       */
      vectorResults.resize(sparsity()->numDerivs());
  

      /* 
       * Carry out the products and sums in the chain rule. 
       *
       * We compute in descending order of derivatives. This lets us
       * reuse the argument's result vector for storage of the 
       * operand's result. 
       *
       */

      /* -- Second derivative terms */

      SUNDANCE_VERB_HIGH(tabs << "evaluating second derivs");

      for (unsigned int i=0; i<d2ArgDerivIndex_.size(); i++)
        {
          const Array<int>& adi = d2ArgDerivIndex_[i];
          const Array<int>& isC = d2ArgDerivIsConstant_[i];

          if (adi[0] == -1)
            {
              vectorResults[d2ResultIndex_[i]] = opDerivs[2]->clone();
              RefCountPtr<EvalVector>& result = vectorResults[d2ResultIndex_[i]];
              if (isC[1]==true)
                {
                  const double& g_u = argConstantResults[adi[1]];
                  if (isC[2]==true)
                    {
                      const double& g_v = argConstantResults[adi[2]];
                      result->multiply_S(g_u * g_v);
                    }
                  else
                    {
                      const EvalVector* g_v = argVectorResults[adi[2]].get();
                      result->multiply_SV(g_u, g_v);
                    }
                }
              else
                {
                  const EvalVector* g_u = argVectorResults[adi[1]].get();
                  if (isC[2]==true)
                    {
                      const double& g_v = argConstantResults[adi[2]];
                      result->multiply_SV(g_v, g_u);
                    }
                  else
                    {
                      const EvalVector* g_v = argVectorResults[adi[2]].get();
                      result->multiply_VV(g_v, g_u);
                    }
                }
            }
          else 
            {
              if (isC[0]==true)
                {
                  vectorResults[d2ResultIndex_[i]] = opDerivs[2]->clone();
                  const double& g_uv = argConstantResults[adi[0]];
                  vectorResults[d2ResultIndex_[i]]->setToConstant(g_uv);
                }
              else
                {
                  vectorResults[d2ResultIndex_[i]] 
                    = argVectorResults[adi[0]];
                }
              const EvalVector* F_g = opDerivs[1].get();
              const EvalVector* F_gg = opDerivs[2].get();
              RefCountPtr<EvalVector>& result = vectorResults[d2ResultIndex_[i]];

              if (isC[1]==true)
                {
                  const double& g_u = argConstantResults[adi[1]];
                  if (isC[2]==true)
                    {
                      const double& g_v = argConstantResults[adi[2]];
                      result->multiply_V_add_SV(F_g, g_u*g_v, F_gg);
                    }
                  else
                    {
                      const EvalVector* g_v = argVectorResults[adi[2]].get();
                      result->multiply_V_add_SVV(F_g, g_u, g_v, F_gg);
                    }
                }
              else
                {
                  const EvalVector* g_u = argVectorResults[adi[1]].get();
                  if (isC[2]==true)
                    {
                      const double& g_v = argConstantResults[adi[2]];
                      result->multiply_V_add_SVV(F_g, g_v, g_u, F_gg);
                    }
                  else
                    {
                      const EvalVector* g_v 
                        = argVectorResults[adi[2]].get();
                      result->multiply_V_add_VVV(F_g, g_v, g_u, F_gg);
                    }
                }
            }
      

        }


      /* --- First derivative terms */
      SUNDANCE_VERB_HIGH(tabs << "evaluating first derivs");

      for (unsigned int i=0; i<d1ArgDerivIndex_.size(); i++)
        {
          if (d1ArgDerivIsConstant_[i])
            {
              const double& g_u = argConstantResults[d1ArgDerivIndex_[i]];
              vectorResults[d1ResultIndex_[i]] = opDerivs[1]->clone();
              vectorResults[d1ResultIndex_[i]]->multiply_S(g_u);
            }
          else
            {
              RefCountPtr<EvalVector>& g_u = argVectorResults[d1ArgDerivIndex_[i]];
              g_u->multiply_V(opDerivs[1].get());
              vectorResults[d1ResultIndex_[i]] = g_u;
            }
        }


      /* --- Zeroth derivative term */

      vectorResults[d0ResultIndex_] = opDerivs[0];
    }  
}



