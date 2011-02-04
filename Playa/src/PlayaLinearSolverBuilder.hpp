/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_LINEARSOLVERBUILDER_HPP
#define PLAYA_LINEARSOLVERBUILDER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;

  /** */
  class LinearSolverBuilder
  {
  public:
    /** */
    static LinearSolver<double> createSolver(const ParameterList& params,
      int verb=0);
    /** */
    static LinearSolver<double> createSolver(const std::string& filename);
  };
  
}

#endif
