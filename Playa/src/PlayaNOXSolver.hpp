/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NOXSOLVER_HPP
#define PLAYA_NOXSOLVER_HPP

#include "PlayaDefs.hpp"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Playa_Group.hpp"
#include "NOX_Playa_StatusTestBuilder.hpp"
#include "NOX_Multiphysics_Solver_Manager.H"
#include "Teuchos_TestForException.hpp"   
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  class NOXSolver 
  {
  public:
    /** */
    NOXSolver(){;}
    /** */
    NOXSolver(const ParameterList& params);
    /** */
    NOXSolver(const ParameterList& nonlinParams,
      const LinearSolver<double>& linSolver);

    /** */
    NOX::StatusTest::StatusType solve(const NonlinearOperator<double>& F, 
      Vector<double>& soln) const ;

    /** */
    const LinearSolver<double>& linSolver() const 
      {return linSolver_;}


  private:

    LinearSolver<double> linSolver_;
    mutable RCP<NOX::StatusTest::Generic> statusTest_;
    mutable ParameterList params_;
    mutable ParameterList printParams_;
  };
}

#endif
