/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_H
#define SUNDANCE_H

/* Utilities */
#include "SundanceDefs.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

/* Symbolics */
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceParameter.hpp"

/* Meshes */
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceMeshTransformation.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceExtrusionMeshTransformation.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceTriangleMeshReader.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"


/* Cell filters */
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceBoundaryCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"

/* Writers */
#include "SundanceFieldWriter.hpp"
#include "SundanceMatlabWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceVTKWriter.hpp"


/* FE  */
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"

/* Problem level classes */
#include "SundanceLinearProblem.hpp"
#include "SundanceL2Projector.hpp"
#include "SundanceNonlinearProblem.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceExprFieldWrapper.hpp"

/* Solvers & stuff */
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFLinearSolverBuilder.hpp"

/* Nonlinear solvers */
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"
#include "TSFNOXSolver.H"


using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceUtils;

namespace SundanceStdFwk
{
  class Sundance
  {
  public:
    static void init(int* argc, void*** argv);
    
    static void finalize();
    
    static void handleException(std::exception& e);


    static void passFailTest(double error, double tol);

    

  private:
    static bool checkTest(double error, double tol);

    static void setSettings(const string& settingsFile);

    static string searchForFile(const string& name);

    static VerbositySetting verbosity(const string& str);

    static string getPathStr();

    static Array<string> parsePathStr();
  };
}




#endif

