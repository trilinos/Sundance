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

#include "SundanceNonlinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace std;
using namespace TSFExtended;


static Time& nlpCtorTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("NonlinearProblem ctor"); 
  return *rtn;
}


NonlinearProblem::NonlinearProblem() 
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(),
    discreteU0_(0)
{
  TimeMonitor timer(nlpCtorTimer());
}


NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
                                   const Expr& eqn, 
                                   const Expr& bc,
                                   const Expr& test, 
                                   const Expr& unk, 
                                   const Expr& u0, 
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    discreteU0_(0)
    
{
  TimeMonitor timer(nlpCtorTimer());

  Expr unkParams;
  Expr fixedParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Expr fixedFieldValues;

  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(test.flattenSpectral()), tuple(unk.flattenSpectral()), tuple(u0),
                          unkParams, unkParamValues,
                          fixedParams, fixedParamValues,
                          tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}

NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
                                   const Expr& eqn, 
                                   const Expr& bc,
                                   const Expr& test, 
                                   const Expr& unk, 
                                   const Expr& u0, 
                                   const Expr& params, 
                                   const Expr& paramVals,  
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    discreteU0_(0)
    
{
  TimeMonitor timer(nlpCtorTimer());

  Expr fixedParams;
  Expr fixedFields;
  Expr fixedParamValues;
  Expr fixedFieldValues;

  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(test.flattenSpectral()), tuple(unk.flattenSpectral()), tuple(u0), 
                          params, paramVals,
                          fixedParams, fixedParamValues,
                          tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}


NonlinearProblem::NonlinearProblem(const RefCountPtr<Assembler>& assembler, 
                                   const Expr& u0)
  : NonlinearOperatorBase<double>(),
    assembler_(assembler),
    u0_(u0),
    discreteU0_(0)
{
  TimeMonitor timer(nlpCtorTimer());
  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem ctor");

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}

TSFExtended::Vector<double> NonlinearProblem::getInitialGuess() const 
{
  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::getInitialGuess()");
  
  Vector<double> u0 = discreteU0_->getVector();

  return u0;
}

void NonlinearProblem::setInitialGuess(const Expr& u0New) 
{
  const DiscreteFunction* in = DiscreteFunction::discFunc(u0New);
  TEST_FOR_EXCEPT(in==0);
  setEvalPt(in->getVector().copy());
  discreteU0_->setVector(currentEvalPt());
}


LinearOperator<double> NonlinearProblem::allocateJacobian() const
{
  return assembler_->allocateMatrix();
}


LinearOperator<double> NonlinearProblem::
computeJacobianAndFunction(Vector<double>& functionValue) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::jacobian()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
                     "null evaluation point in "
                     "NonlinearProblem::jacobian()");

  discreteU0_->setVector(currentEvalPt());

  assembler_->assemble(J_, functionValue);

  return J_;
}

void NonlinearProblem::
computeJacobianAndFunction(LinearOperator<double>& J,
                           Vector<double>& resid) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::jacobian()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
                     "null evaluation point in "
                     "NonlinearProblem::jacobian()");

  TEST_FOR_EXCEPTION(J.ptr().get()==0, RuntimeError,
                     "null Jacobian pointer in "
                     "NonlinearProblem::jacobian()");

  TEST_FOR_EXCEPTION(resid.ptr().get()==0, RuntimeError,
                     "null residual pointer in "
                     "NonlinearProblem::jacobian()");

  discreteU0_->setVector(currentEvalPt());

  assembler_->assemble(J, resid);
  J_ = J;
}


Vector<double> NonlinearProblem::computeFunctionValue() const 
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::computeFunctionValue()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
                     "null evaluation point in "
                     "NonlinearProblem::computeFunctionValue()");

  discreteU0_->setVector(currentEvalPt());

  Vector<double> rtn = createMember(range());
  assembler_->assemble(rtn);

  return rtn;
}



void NonlinearProblem::computeFunctionValue(Vector<double>& resid) const 
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::computeFunctionValue()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
                     "null evaluation point in "
                     "NonlinearProblem::computeFunctionValue()");

  TEST_FOR_EXCEPTION(resid.ptr().get()==0, RuntimeError,
                     "null residual pointer in "
                     "NonlinearProblem::computeFunctionValue()");

  discreteU0_->setVector(currentEvalPt());

  assembler_->assemble(resid);
}

