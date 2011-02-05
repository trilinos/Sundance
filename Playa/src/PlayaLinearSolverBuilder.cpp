/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaExceptions.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaAmesosSolver.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaBelosSolver.hpp"
#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaBlockTriangularSolverDecl.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaBICGSTABSolverImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"
#endif

using namespace Playa;
using namespace PlayaOps;
using namespace Teuchos;


LinearSolver<double> LinearSolverBuilder::createSolver(const std::string& filename)
{
  ParameterXMLFileReader reader(filename);
  ParameterList solverParams = reader.getParameters();
  return createSolver(solverParams);
}



LinearSolver<double> LinearSolverBuilder::createSolver(const ParameterList& params, int verb)
{
  TEST_FOR_EXCEPTION(!params.isSublist("Linear Solver"), std::runtime_error,
    "did not find Linear Solver sublist in " << params);

  ParameterList solverSublist = params.sublist("Linear Solver");

  const std::string& solverType = getParameter<string>(solverSublist, "Type");

  Tabs tab;
  PLAYA_MSG1(verb, tab << "Solver builder creating a solver of type="
    << solverType);
  Tabs tab2;
  PLAYA_MSG2(verb, tab2 << "params = " << solverSublist);

  if (solverType=="Aztec")
  {
    return new AztecSolver(solverSublist);
  }
  else if (solverType=="Playa" || solverType=="TSF")
  {
    const std::string& solverMethod = getParameter<string>(solverSublist, "Method");
    if (solverMethod=="BICGSTAB") 
    {
      return new BICGSTABSolver<double>(solverSublist);
    }
    else if (solverMethod=="GMRES")
    {
      TEST_FOR_EXCEPTION(true, RuntimeError, "Playa GMRES solver not implemented");
    }
  }
  else if (solverType=="Amesos")
  {
    return new AmesosSolver(solverSublist);
  }
  else if (solverType=="Belos")
  {
    return new BelosSolver(solverSublist);
  }
  else if (solverType=="Block Triangular")
  {
    ParameterList subSolverParams = solverSublist.sublist("Sub Solver");
    LinearSolver<double> subSolver = createSolver(subSolverParams);
    return new BlockTriangularSolver<double>(subSolver);
  }

  TEST_FOR_EXCEPTION(true, std::runtime_error, 
    "Could not create a solver from parameter list " 
    << params);
  return LinearSolver<double>();
    
}

