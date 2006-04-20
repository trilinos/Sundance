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




#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceNLPModelEvaluator.hpp"

#ifdef HAVE_ENABLED_MOOCHO

SundanceNLPModelEvaluator::SundanceNLPModelEvaluator(const VectorType<double>& vecType)
  : SundanceModelEvaluator(vecType),
    paramSpace_(),
    paramExpr_(),
    stateExpr_(),
    prob_(),
    sensProb_(),
    obj_(),
    objEval_()
{;}


void SundanceNLPModelEvaluator::initialize(const Expr& paramExpr,
                                           const Expr& stateExpr,
                                           const Expr& stateExprVal,
                                           const NonlinearProblem& prob, 
                                           const Array<LinearProblem>& sensProb,
                                           const Functional& objective)
{
  paramExpr_ = paramExpr.flatten();
  stateExpr_ = stateExprVal;
  prob_ = prob;
  sensProb_ = sensProb;
  obj_ = objective;
  objEval_ = obj_.evaluator(stateExpr, stateExprVal);
  RefCountPtr<const Thyra::VectorSpaceBase<double> > pSpace
    = rcp(new Thyra::DefaultSerialVectorSpace<double>(paramExpr_.size()));
  paramSpace_ = pSpace;
}

void SundanceNLPModelEvaluator::internalEvalModel(const Vector<double>& stateVec,
                                                  const Vector<double>& params,
                                                  Vector<double>& resid,
                                                  double& objFuncVal,
                                                  LinearOperator<double>& df_dx,
                                                  Array<Vector<double> >& df_dp,
                                                  Vector<double>& dg_dp_T,
                                                  Vector<double>& dg_dx_T) const 
{
  Tabs tabs;
  TEST_FOR_EXCEPTION(params.ptr().get()==0, RuntimeError,
                     "null params vector!");
  TEST_FOR_EXCEPTION( (int) paramExpr_.size() != params.space().dim(),
                      RuntimeError,
                      "Mismatch between input parameter vector space and "
                      "symbolic parameter object size");

  TEST_FOR_EXCEPTION(stateVec.ptr().get()==0, RuntimeError,
                     "null state vector!");


  /* set the symbolic parameter values to the input parameters */
  if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "parameters are: ")
    for (unsigned int i=0; i<paramExpr_.size(); i++)
      {
        Tabs tab2;
        Expr p_i = paramExpr_[i];
        p_i.setParameterValue(params[i]);
        if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << i << " " << paramExpr_[i]);
      }

      

  /* set the symbolic state function value to the input state */

  if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "state vector is " 
                       << stateVec.description());
  prob_.setEvalPt(stateVec);
  

  /* compute the Jacobian and residual of the constraints */
  if (df_dx.ptr().get()!=0 && resid.ptr().get() != 0)
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing Jacobian and residual...");
      prob_.computeJacobianAndFunction(df_dx, resid);          
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_HIGH(tabs << "residual vector is " 
                         << resid.description());
    }
  else if (resid.ptr().get() != 0 && df_dx.ptr().get()==0)
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing residual w/o Jacobian...");
      prob_.computeFunctionValue(resid);         
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_HIGH(tabs << "residual vector is " 
                         << resid.description());
    }
  else if (df_dx.ptr().get()!=0 && resid.ptr().get() == 0)
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing Jacobian w/o residual...");
      Vector<double> dummy = constraintSpace().createMember();
      prob_.computeJacobianAndFunction(df_dx, dummy);          
    }
  else
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "requested neither the residual "
                           "nor the Jacobian");
    }
        

  /* compute the derivatives of the constraint residual wrt the parameters */
  if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing sensitivities...");
  for (unsigned int i=0; i<paramExpr_.size(); i++)
    {
      df_dp[i] = sensProb_[i].getRHS();
    }

  /* evaluate the objective function and its gradient wrt the states */ 
  if (dg_dx_T.ptr().get() != 0)
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing objective function "
                           "and gradient");
      Expr df_dx_Expr = objEval_.evalGradient(objFuncVal);
      dg_dx_T.acceptCopyOf(DiscreteFunction::discFunc(df_dx_Expr)->getVector());
    }
  else
    {
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "computing objective function");
      objFuncVal = objEval_.evaluate();
    }
      


  /* evaluate the derivatives of the objective function wrt 
   * the parameters */
  /* set the symbolic parameter values to the input parameters */
  if (dg_dp_T.ptr().get() != 0)
    {
      dg_dp_T.setToConstant(0.0);
    }
}


Array<double> SundanceNLPModelEvaluator::parameters() const
{
  Array<double> rtn;
  for (unsigned int i=0; i<paramExpr_.size(); i++)
    {
      const SpatiallyConstantExpr* sce 
        = dynamic_cast<const SpatiallyConstantExpr*>(paramExpr_[i].ptr().get());
      rtn.append(sce->value());
    }
  return rtn;
}



#endif

