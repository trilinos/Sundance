#ifndef PDEOPT_LINEARPDECONSTRAINEDOBJ_H
#define PDEOPT_LINEARPDECONSTRAINEDOBJ_H

#include "PDEOptPDEConstrainedObjBase.hpp"
#include "SundanceExpr.hpp"
#include "SundanceFunctional.hpp"
#include "SundanceLinearProblem.hpp"

namespace Sundance
{
using namespace Playa;


/**
 * LinearPDEConstrainedObj is a base class for objective functions of the
 * reduced-space variable where the constraint is a linear PDE in the
 * state variables. 
 */
class LinearPDEConstrainedObj : public PDEConstrainedObjBase
{
public:
  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const LinearSolver<double>& solver,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const Array<LinearSolver<double> >& solvers,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const LinearSolver<double>& solver,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const Array<LinearSolver<double> >& solvers,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /** virtual dtor */
  virtual ~LinearPDEConstrainedObj(){;}



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

  Array<LinearProblem> stateProbs_;

  Array<LinearProblem> adjointProbs_;

  Array<LinearSolver<double> > solvers_;

};

}

#endif
