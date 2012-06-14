/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaExceptions.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaNOXSolver.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif

using namespace Playa;
using namespace PlayaExprTemplates;
using namespace Teuchos;


NonlinearSolver<double> NonlinearSolverBuilder::createSolver(const std::string& filename)
{
  ParameterXMLFileReader reader(filename);
  ParameterList solverParams = reader.getParameters();
  return createSolver(solverParams);
}



NonlinearSolver<double> NonlinearSolverBuilder::createSolver(const ParameterList& params)
{
  if (params.isSublist("NOX Solver"))
  {
    return new NOXSolver(params);
  }
  else if (params.isSublist("Nonlinear Solver"))
  {
    ParameterList sub = params.sublist("Nonlinear Solver");
    Array<string> names = tuple<string>("Newton Armijo Solver", "Newton-Armijo Solver", "NewtonArmijoSolver");
    for (int i=0; i<names.size(); i++)
    {
      if (sub.isSublist(names[i]))
      {
        ParameterList subsub = sub.sublist(names[i]);
        LinearSolver<double> linSolver;
        if (subsub.isParameter("Linear Solver"))
        {
          string solverFile = subsub.get<string>("Linear Solver");
          linSolver = LinearSolverBuilder::createSolver(solverFile);
        }
        else if (subsub.isSublist("Linear Solver"))
        {
          linSolver = LinearSolverBuilder::createSolver(subsub);
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
            "Nonlinear solver parameter list " << sub
            << " does not appear to specify a solver for the linear subproblems");
        }
        return new NewtonArmijoSolver<double>(subsub, linSolver);
      }
    }
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
      "Nonlinear solver parameter list " << params
      << " can't be parsed to find a nonlinear solver");
  }

  return NonlinearSolver<double>();
    
}

