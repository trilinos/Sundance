#ifndef PLAYA_PCG_SOLVER_HPP
#define PLAYA_PCG_SOLVER_HPP

#include "PlayaLinearSolverBaseImpl.hpp"
#include "PlayaPreconditionerFactory.hpp"

namespace Playa
{

  /** */
  class PCGSolver : public LinearSolverBase<double>,
			 public Handleable<LinearSolverBase<double> >
  {
  public:
    PCGSolver(const ParameterList& params,
	      const PreconditionerFactory<double>& precFactory)
      : LinearSolverBase<double>(params), precFactory_(precFactory) {}

    SolverState<double> solve(const LinearOperator<double>& A,
			      const Vector<double>& b,
			      Vector<double>& x) const ;

    SolverState<double> solveUnprec(const LinearOperator<double>& A,
				    const Vector<double>& b,
				    Vector<double>& x) const ;

    RCP<LinearSolverBase<double> > getRcp() {return rcp(this);}

  private:
    PreconditionerFactory<double> precFactory_;
  };
}


#endif
