/* @HEADER@ */
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

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;



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
    rhs_()
{
  Expr u = unk.flatten();
  Expr v = test.flatten();
  cerr << "unk  = " << u << endl;
  cerr << "unk size = " << u.size() << endl;
  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
    {
      Expr z = new ZeroExpr();
      zero[i] = z;
    }

  Expr u0 = new ListExpr(zero);
  cerr << "zero = " << u0 << endl;
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0));

  assembler_ = rcp(new Assembler(mesh, eqnSet, vecType));
}

TSFExtended::Vector<double> LinearProblem::getRHS() const 
{
  assembler_->assemble(A_, rhs_);
  return rhs_;
}


TSFExtended::LinearOperator<double> LinearProblem::getOperator() const 
{
  assembler_->assemble(A_, rhs_);
  return A_;
}

Expr LinearProblem::solve(const LinearSolver<double>& solver) const 
{
  Vector<double> solnVec;

  assembler_->assemble(A_, rhs_);
  rhs_.scale(-1.0);

  SolverState<double> state = solver.solve(A_, rhs_, solnVec);
  TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
                     RuntimeError,
                     "solve failed!");

  

  Expr soln = new DiscreteFunction(*(assembler_->solutionSpace()),
                                   solnVec, "soln");

  return soln;
}

