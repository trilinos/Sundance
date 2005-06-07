/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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
#include "Teuchos_CommandLineProcessor.hpp"

/* Symbolics */
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceParameter.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"

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
#include "SundanceMeshBuilder.hpp"
#include "SundanceBamgMeshReader.hpp"


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
#include "SundanceFunctional.hpp"
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

/* do explicit qualification of List to avoid conflicts
 * with the unfriendly lack of namespaces in MPI2C++
 */
using SundanceCore::List;



namespace SundanceStdFwk
{
  /**
   * Class Sundance provides several static methods for
   * managing the environment of a simulation. Every simulation code
   * should begin with a call to Sundance::init() and end with
   * a call to Sundance::finalize().
   */
  class Sundance
  {
  public:
    static void setOption(const string& optionName, 
                          int& value, 
                          const string& helpMsg);

    static void setOption(const string& optionName, 
                          string& value, 
                          const string& helpMsg);

    static void setOption(const string& optionName, 
                          double& value, 
                          const string& helpMsg);

    static void setOption(const string& optionTrueName, 
                          const string& optionFalseName, 
                          bool& value, 
                          const string& helpMsg);

    /** 
     * Do initialization steps such as starting MPI (if necessary), 
     * parsing the Unix command
     * line, and reading runtime options from the configuration file.
     * MPI is initialized through a call to Teuchos::MPISession::init(), 
     * which in turn checks whether MPI needs initialization and calls
     * MPI_Init() if necessary. If some other library you're using has
     * its own MPI initialization system, it is thus perfectly safe to
     * call Sundance::init() as well.
     */
    static void init(int* argc, void*** argv);
    
    /** 
     * Do finalization steps such as calling MPI_finalize() (if necessary),
     * and collecting and printing timing information.
     */
    static void finalize();
    
    /** */
    static void handleException(std::exception& e);

    /** */
    static void passFailTest(double error, double tol);

    /** */
    static string searchForFile(const string& name);    

    static VerbositySetting verbosity(const string& str);

    /** */
    static void setSettings(const XMLObject& xml);

  private:
    static CommandLineProcessor& clp()
    {static CommandLineProcessor rtn; return rtn;}

    static bool checkTest(double error, double tol);

    static void setSettings(const string& settingsFile);




    static string getPathStr();

    static Array<string> parsePathStr();
  };
}




#endif

