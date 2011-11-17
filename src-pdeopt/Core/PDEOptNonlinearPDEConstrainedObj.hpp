#ifndef PDEOPT_NONLINEARPDECONSTRAINEDOBJ_H
#define PDEOPT_NONLINEARPDECONSTRAINEDOBJ_H

#include "PDEOptPDEConstrainedObjBase.hpp"
#include "SundanceExpr.hpp"
#include "PlayaNOXSolver.hpp"
#include "SundanceFunctional.hpp"
#include "SundanceNonlinearProblem.hpp"

namespace Sundance
{
using namespace Playa;


/**
 * NonlinearPDEConstrainedObj is a base class for objective functions of the
 * reduced-space variable where the constraint is a nonlinear PDE in the
 * state variables. 
 */
class NonlinearPDEConstrainedObj : public PDEConstrainedObjBase
{
public:
  /**  */
  NonlinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const NOXSolver& solver,
    const LinearSolver<double>& adjSolver,
    int verb=0);

  /**  */
  NonlinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const NOXSolver& solver,
    const LinearSolver<double>& adjSolver,
    int verb=0);

  /**  */
  NonlinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const NOXSolver& solver,
    const LinearSolver<double>& adjSolver,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /**  */
  NonlinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const NOXSolver& solver,
    const LinearSolver<double>& adjSolver,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /** virtual dtor */
  virtual ~NonlinearPDEConstrainedObj(){;}



  /** Solve the sequence of state equations, followed by postprocessing.
   * At the end of this call, the system is ready for evaluation of
   * the objective function or solution of the adjoint equations. */
  void solveState(const Vector<double>& x) const;

  /** Solve the sequence of state equations, then do postprocessing,
   * then finally the adjoint equations in <b>reverse</b> order. 
   * At the end of this call, the system is ready for evaluation
   * of the objective function and its gradient. */
  void solveStateAndAdjoint(const Vector<double>& x) const;

  /** Set up the linear equations */
  void initEquations(
    const Array<Expr>& stateVars,
    const Array<Expr>& adjointVars,
    const Array<Array<Expr> >& fixedVarsInStateEqns,
    const Array<Array<Expr> >& fixedVarsInStateEqnsVals,
    const Array<Array<Expr> >& fixedVarsInAdjointEqns,
    const Array<Array<Expr> >& fixedVarsInAdjointEqnsVals
    );


private:

  Array<NonlinearProblem> stateProbs_;

  Array<LinearProblem> adjointProbs_;

  NOXSolver solver_;

  LinearSolver<double> adjSolver_;
};

}

#endif
