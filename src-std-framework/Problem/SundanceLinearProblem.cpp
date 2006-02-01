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
  for (unsigned int i=0; i<u.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr u0 = new ListExpr(zero);
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(v), tuple(u), tuple(u0)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType)));
}

LinearProblem::LinearProblem(const Mesh& mesh, 
                             const Expr& eqn, 
                             const Expr& bc,
                             const BlockArray& test, 
                             const BlockArray& unk)
  : assembler_(),
    A_(),
    rhs_(),
    status_()
{
  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0(unk.size());
  Array<VectorType<double> > testVecType(test.size());
  Array<VectorType<double> > unkVecType(unk.size());

  for (unsigned int i=0; i<test.size(); i++)
    {
      v[i] = test[i].expr().flatten();
      testVecType[i] = test[i].vecType();
    }

  for (unsigned int i=0; i<unk.size(); i++)
    {
      u[i] = unk[i].expr().flatten();
      Array<Expr> zero(u[i].size());
      for (unsigned int j=0; j<u[i].size(); j++) 
        {
          Expr z = new ZeroExpr();
          zero[j] = z;
        }
      u0[i] = new ListExpr(zero);
    }
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType));
}

LinearProblem::LinearProblem(const RefCountPtr<Assembler>& assembler) 
  : assembler_(assembler),
    A_(),
    rhs_()
{}


TSFExtended::Vector<double> LinearProblem::getRHS() const 
{
  Tabs tab;
  SUNDANCE_VERB_MEDIUM(tab << "LinearProblem::solve() building vector");
  assembler_->assemble(rhs_);
  return rhs_;
}


TSFExtended::LinearOperator<double> LinearProblem::getOperator() const 
{
  Tabs tab;
  SUNDANCE_VERB_MEDIUM(tab << "LinearProblem::solve() building matrix and vector");
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
  
  SUNDANCE_VERB_MEDIUM(tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  rhs_.scale(-1.0);

  SUNDANCE_VERB_MEDIUM(tab << "LinearProblem::solve() solving system");

  status_ = rcp(new SolverState<double>(solver.solve(A_, rhs_, solnVec)));

  const SolverState<double>& state = *status_;
  SUNDANCE_VERB_MEDIUM(tab << 
                       "LinearProblem::solve() done solving system: status is " 
                       << state.stateDescription());

  Expr soln;
  if (state.finalState() != SolveConverged) 
    {
      handleSolveFailure();
    }
  else
    {
      soln = formSolutionExpr(solnVec);
    }

  return soln;
}

SolverState<double> LinearProblem
::solve(const LinearSolver<double>& solver,
        Expr& soln) const 
{
  Tabs tab;
  Vector<double> solnVec;
  
  SUNDANCE_VERB_MEDIUM(tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  rhs_.scale(-1.0);

  SUNDANCE_VERB_LOW(tab << "solving LinearProblem");

  status_ = rcp(new SolverState<double>(solver.solve(A_, rhs_, solnVec)));

  const SolverState<double>& state = *status_;
  SUNDANCE_VERB_MEDIUM(tab << 
                       "LinearProblem::solve() done solving system: status is " 
                       << state.stateDescription());

  if (state.finalState() != SolveConverged) 
    {
      handleSolveFailure();
    }
  else
    {
      if (soln.ptr().get()==0)
        {
          soln = formSolutionExpr(solnVec);
        }
      else
        {
          DiscreteFunction::discFunc(soln)->setVector(solnVec);
        }
    }

  return state;
}

Expr LinearProblem::formSolutionExpr(const Vector<double>& solnVector) const
{
  Array<Expr> rtn(assembler_->solutionSpace().size());
  for (int i=0; i<rtn.size(); i++)
    {
      string name = "soln";
      if ((int)rtn.size() > 1) name += "[" + Teuchos::toString(i) + "]";
      rtn[i] = new DiscreteFunction(*(assembler_->solutionSpace()[i]),
                                    solnVector.getBlock(i), name);
    }
  if ((int) rtn.size() > 1)
    {
      return new ListExpr(rtn);
    }
  else
    {
      return rtn[0];
    }
}

void LinearProblem::handleSolveFailure() const 
{
  const SolverState<double>& state = *status_;

  TeuchosOStringStream ss;
  ss << "Solve failed! state = "
     << state.stateDescription()
     << "\nmessage=" << state.finalMsg()
     << "\niters taken = " << state.finalIters()
     << "\nfinal residual = " << state.finalResid();
  
  if (dumpBadMatrix())
    {
      if (A_.ptr().get() != 0)
        {
          ofstream osA(badMatrixFilename().c_str());
          A_.print(osA);
          ss << "\nmatrix written to " << badMatrixFilename();
        }
      else
        {
          ss << "\nthe matrix is null! Evil is afoot in your code...";
        }
      if (rhs_.ptr().get() != 0)
        {
          ofstream osb(badVectorFilename().c_str());
          rhs_.print(osb);
          ss << "\nRHS vector written to " << badVectorFilename();
        }
      else
        {
          ss << "\nthe RHS vector is null! Evil is afoot in your code...";
        }
    }

  TEST_FOR_EXCEPTION((state.finalState() != SolveConverged) && stopOnSolveFailure(),
                     RuntimeError, TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss));

  cerr << TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss) << endl;
}
