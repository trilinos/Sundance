#include "PDEOptLinearPDEConstrainedObj.hpp"
#include "Sundance.hpp"

using namespace Teuchos;
using namespace Playa;

namespace Sundance
{

LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Expr& stateVars,
  const Expr& stateVarVals,
  const Expr& adjointVars,
  const Expr& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const LinearSolver<double>& solver,
  const RCP<IterCallbackBase>& iterCallback,
  int verb)
  : PDEConstrainedObjBase(lagrangian, tuple(stateVarVals),
    tuple(adjointVarVals), designVarVal, iterCallback, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(tuple(solver))
{
  init(tuple(stateVars), tuple(adjointVars), designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Array<Expr>& stateVars,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVars,
  const Array<Expr>& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const Array<LinearSolver<double> >& solvers,
  const RCP<IterCallbackBase>& iterCallback,
  int verb)
  : PDEConstrainedObjBase(lagrangian, stateVarVals,
    adjointVarVals, designVarVal, iterCallback, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(solvers)
{
  init(stateVars, adjointVars, designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Expr& stateVars,
  const Expr& stateVarVals,
  const Expr& adjointVars,
  const Expr& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const LinearSolver<double>& solver,
  int verb)
  : PDEConstrainedObjBase(lagrangian, tuple(stateVarVals),
    tuple(adjointVarVals), designVarVal, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(tuple(solver))
{
  init(tuple(stateVars), tuple(adjointVars), designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Array<Expr>& stateVars,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVars,
  const Array<Expr>& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const Array<LinearSolver<double> >& solvers,
  int verb)
  : PDEConstrainedObjBase(lagrangian, stateVarVals,
    adjointVarVals, designVarVal, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(solvers)
{
  init(stateVars, adjointVars, designVar);
}


void LinearPDEConstrainedObj::initEquations(
  const Array<Expr>& stateVars,
  const Array<Expr>& adjointVars,
  const Array<Array<Expr> >& fixedVarsInStateEqns,
  const Array<Array<Expr> >& fixedVarsInStateEqnsVals,
  const Array<Array<Expr> >& fixedVarsInAdjointEqns,
  const Array<Array<Expr> >& fixedVarsInAdjointEqnsVals
  )
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "setting up linear equations");
  
  for (int i=0; i<stateVars.size(); i++)
  {
    Tabs tab1;
    PLAYA_MSG3(verb(), tab1 << "setting up state equation #" << i);
    Expr fixedVars = new ListExpr(fixedVarsInStateEqns[i]);
    Expr fixedVarVals = new ListExpr(fixedVarsInStateEqnsVals[i]);
    LinearProblem stateProb 
      = Lagrangian().linearVariationalProb(adjointVars[i], adjointVarVals(i),
        stateVars[i],
        fixedVars, fixedVarVals);
                                   
    stateProbs_.append(stateProb);
  }

  for (int i=0; i<adjointVars.size(); i++)
  {
    Tabs tab1;
    PLAYA_MSG3(verb(), tab1 << "setting up adjoint equation #" << i);
    Expr fixedVars = new ListExpr(fixedVarsInAdjointEqns[i]);
    Expr fixedVarVals = new ListExpr(fixedVarsInAdjointEqnsVals[i]);
    LinearProblem adjointProb 
      = Lagrangian().linearVariationalProb(stateVars[i], stateVarVals(i),
        adjointVars[i],
        fixedVars, fixedVarVals);
                                   
    adjointProbs_.append(adjointProb);
  }

  PLAYA_MSG2(verb(), tab << "done setting up linear equations");
}




void LinearPDEConstrainedObj::solveState(const Vector<double>& x) const
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "solving state"); 
  PLAYA_MSG3(verb(), tab << "|x|=" << x.norm2()); 
  PLAYA_MSG5(verb(), tab << "x=" << endl << tab << x.norm2());
  setDiscreteFunctionVector(designVarVal(), x);

  /* solve the state equations in order */
  for (int i=0; i<stateProbs_.size(); i++)
  {
    SolverState<double> status 
      = stateProbs_[i].solve(solvers_[i], stateVarVals(i));
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "state equation could not be solved: status="
      << status.stateDescription());
  }

  PLAYA_MSG2(verb(), tab << "done state solve"); 
  /* do postprocessing */
  statePostprocCallback();
}



void LinearPDEConstrainedObj
::solveStateAndAdjoint(const Vector<double>& x) const
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "solving state and adjoint"); 
  PLAYA_MSG3(verb(), tab << "|x|=" << x.norm2()); 
  PLAYA_MSG5(verb(), tab << "x=" << endl << tab << x.norm2()); 

  Tabs tab1;
  setDiscreteFunctionVector(designVarVal(), x);

  PLAYA_MSG3(verb(), tab1 << "solving state eqns");
  /* solve the state equations in order */
  for (int i=0; i<stateProbs_.size(); i++)
  {
    SolverState<double> status 
      = stateProbs_[i].solve(solvers_[i], stateVarVals(i));

    /* if the solve failed, write out the design var and known state
     * variables */
    if (status.finalState() != SolveConverged)
    {
      FieldWriter w = new VTKWriter("badSolve");
      w.addMesh(Lagrangian().mesh());
      w.addField("designVar", new ExprFieldWrapper(designVarVal()));
      for (int j=0; j<i; j++)
      {
        Expr tmp = stateVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("stateVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }
      }
      w.write();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "state equation " << i 
      << " could not be solved: status="
      << status.stateDescription());
  }

  PLAYA_MSG3(verb(), tab1 << "done solving state eqns");

  /* do postprocessing */
  statePostprocCallback();

  PLAYA_MSG3(verb(), tab1 << "solving adjoint eqns");

  /* solve the adjoint equations in reverse order */
  for (int i=adjointProbs_.size()-1; i>=0; i--)
  {
    SolverState<double> status 
      = adjointProbs_[i].solve(solvers_[i], adjointVarVals(i));

    /* if the solve failed, write out the design var and known state
     * and adjoint variables */
    if (status.finalState() != SolveConverged)
    {
      FieldWriter w = new VTKWriter("badSolve");
      w.addMesh(Lagrangian().mesh());
      w.addField("designVar", new ExprFieldWrapper(designVarVal()));
      for (int j=0; j<stateProbs_.size(); j++)
      {
        Expr tmp = stateVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("stateVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }
      }
      for (int j=adjointProbs_.size()-1; j>i; j--)
      {
        Expr tmp = adjointVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("adjointVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }

      }
      w.write();

    }
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "adjoint equation " << i 
      << " could not be solved: status="
      << status.stateDescription());
  }
  PLAYA_MSG3(verb(), tab1 << "done solving adjoint eqns");
  PLAYA_MSG2(verb(), tab1 << "done solving state and adjoint eqns");
}

  


}










