/* @HEADER@ */
/* @HEADER@ */

#include "SundanceLinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceZeroExpr.hpp"

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
  Expr u0 = new ZeroExpr();
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, test, unk, u0, 
                          rcp(new BruteForceEvaluatorFactory())));

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
  
  Expr soln = new DiscreteFunction(*(assembler_->solutionSpace()),
                                   solnVec, "u0");

  return soln;
}

