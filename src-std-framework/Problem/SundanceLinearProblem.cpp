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

#include "SundanceLinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceListExpr.hpp"
#include "TSFSolverState.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;
using namespace std;


LinearProblem::LinearProblem() 
  : assembler_(),
    A_(),
    rhs_()
{}


LinearProblem::LinearProblem(const Mesh& mesh, 
                                   const Expr& eqn, 
                                   const Expr& bc,
                                   const Expr& test, 
                                   const Expr& unk, 
                                   const VectorType<double>& vecType)
  : assembler_(),
    A_(),
    rhs_(),
    status_()
{
  Expr u = unk.flatten();
  Expr v = test.flatten();
  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr u0 = new ListExpr(zero);
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0));

  assembler_ = rcp(new Assembler(mesh, eqnSet, vecType));
}

LinearProblem::LinearProblem(const RefCountPtr<Assembler>& assembler) 
  : assembler_(assembler),
    A_(),
    rhs_()
{}


TSFExtended::Vector<double> LinearProblem::getRHS() const 
{
  Tabs tab;
  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() building vector");
  assembler_->assemble(rhs_);
  return rhs_;
}


TSFExtended::LinearOperator<double> LinearProblem::getOperator() const 
{
  Tabs tab;
  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() building matrix and vector");
  assembler_->assemble(A_, rhs_);
  return A_;
}

SolverState<double> LinearProblem::solveStatus() const
{
  TEST_FOR_EXCEPTION(status_.get()==0, RuntimeError,
                     "Null status pointer in LinearProblem::solveStatus().");
  return *status_;
}

Expr LinearProblem::solve(const LinearSolver<double>& solver) const 
{
  Tabs tab;
  Vector<double> solnVec;
  
  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  rhs_.scale(-1.0);

  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() solving system");

  status_ = rcp(new SolverState<double>(solver.solve(A_, rhs_, solnVec)));

  const SolverState<double>& state = *status_;
  SUNDANCE_VERB_LOW(tab << 
                    "LinearProblem::solve() done solving system: status is " 
                    << state.stateDescription());

  Expr soln = formSolutionExpr(solnVec);

  return soln;
}

SolverState<double> LinearProblem
::solve(const LinearSolver<double>& solver,
        Expr& soln) const 
{
  Tabs tab;
  Vector<double> solnVec;
  
  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  rhs_.scale(-1.0);

  SUNDANCE_VERB_LOW(tab << "LinearProblem::solve() solving system");

  status_ = rcp(new SolverState<double>(solver.solve(A_, rhs_, solnVec)));

  const SolverState<double>& state = *status_;
  SUNDANCE_VERB_LOW(tab << 
                    "LinearProblem::solve() done solving system: status is " 
                    << state.stateDescription());

  if (soln.ptr().get()==0)
    {
      soln = formSolutionExpr(solnVec);
    }
  else
    {
      DiscreteFunction::discFunc(soln)->setVector(solnVec);
    }

  return state;
}

Expr LinearProblem::formSolutionExpr(const Vector<double>& solnVector) const
{
  return new DiscreteFunction(*(assembler_->solutionSpace()),
                              solnVector, "soln");
}

