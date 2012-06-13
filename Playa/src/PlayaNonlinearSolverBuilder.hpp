/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NONLINEARSOLVERBUILDER_HPP
#define PLAYA_NONLINEARSOLVERBUILDER_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;

  /** */
  class NonlinearSolverBuilder
  {
  public:
    /** */
    static NonlinearSolver<double> createSolver(const ParameterList& params);
    /** */
    static NonlinearSolver<double> createSolver(const std::string& filename);
  };
  
}

#endif
