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

#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFLinearCombination.hpp"

using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;});

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

        int np = MPIComm::world().getNProc();

      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);

      mesher.serializeLocal();

      Mesh mesh = mesher.getMesh();

      Array<int> funcs = tuple(0);

      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      
      Expr x = new CoordExpr(0);
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      //      Expr u0 = new DiscreteFunction(new Lagrange(1), "u0");
      Expr u0 = new ZeroExpr();
      Expr v = new TestFunction(new Lagrange(2), "v");
      Expr dx = new Derivative(0);

      Expr exactSoln = -x*x + 2.0*x;

      QuadratureFamily quad = new GaussianQuadrature(2);
      
      Expr eqn = Integral(interior, (dx*v)*(dx*u) - 2.0*v, quad);
      Expr bc = EssentialBC(leftPoint, v*u, quad);

      RefCountPtr<EquationSet> eqnSet 
        = rcp(new EquationSet(eqn, bc, v, u, u0));

      //      verbosity<EvalVector>() = VerbExtreme;
      // verbosity<QuadratureEvalMediator>() = VerbExtreme;
      //verbosity<EvaluatableExpr>() = VerbExtreme;
      //verbosity<Evaluator>() = VerbExtreme;
      //      EvalVector::shadowOps() = true;

      //Evaluator::classVerbosity() = VerbHigh;
//       Assembler::classVerbosity() = VerbExtreme;
//       IntegralGroup::classVerbosity() = VerbHigh;
//       Expr::showAllParens() = true;

      VectorType<double> vecType = new EpetraVectorType();

      Assembler assembler(mesh, eqnSet, vecType); 

      LinearOperator<double> A;
      Vector<double> b;
      Vector<double> solnVec;

      assembler.assemble(A, b);

      ParameterList solverParams;

      solverParams.set(LinearSolverBase<double>::verbosityParam(), 2);
      solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
      solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-12);

      LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);

      SolverState<double> state = solver.solve(A, b, solnVec);

      cerr << "solver state = " << endl << state << endl;

      Expr soln = new DiscreteFunction(*(assembler.solutionSpace()),
                                       solnVec, "u0");

      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              new GaussianQuadrature(6));

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
                                   derivErr*derivErr, 
                                   new GaussianQuadrature(4));

      Expr junk = dx*exactSoln - soln[0];
      
      Expr junkExpr = Integral(interior, 
                               junk*junk,
                               new GaussianQuadrature(12));

      FunctionalEvaluator errInt(mesh, errExpr);
      FunctionalEvaluator junkInt(mesh, junkExpr);
      FunctionalEvaluator derivErrInt(mesh, derivErrExpr);
      FunctionalEvaluator exactInt(mesh, Integral(interior,exactSoln,quad));
      FunctionalEvaluator solnInt(mesh, Integral(interior,soln[0], quad));

      //EvaluatableExpr::classVerbosity() = VerbHigh;
      FunctionalEvaluator dExactInt(mesh, Integral(interior,dx*exactSoln,
                                                   new GaussianQuadrature(8)));

      cerr << endl << endl;
      cerr << "###################################################"
           << endl << endl;

      FunctionalEvaluator dSolnInt(mesh, Integral(interior,dx*soln[0], 
                                                  new GaussianQuadrature(10)));
      EvaluatableExpr::classVerbosity() = VerbSilent;
      double errorSq = errInt.evaluate();
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double derivErrorSq = derivErrInt.evaluate();
      cerr << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;

      cerr << "integral(exact soln) = " << exactInt.evaluate() << endl << endl;
      cerr << "integral(numerical soln) = " << solnInt.evaluate() << endl << endl;

      cerr << "integral(dx*exact soln) = " << dExactInt.evaluate() << endl << endl;

      //Evaluator::classVerbosity() = VerbHigh;
      cerr << "integral(dx*numerical soln) = " << dSolnInt.evaluate() << endl << endl;
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
