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

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceExodusMeshReader.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;}) 
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 

#if defined(HAVE_SUNDANCE_EXODUS)

int main(int argc, char** argv)
{
  
  try
		{
      int nx = 8;
      int ny = 8;
      string meshFile="builtin";
      string solverFile = "aztec-ml.xml";
      Sundance::setOption("meshFile", meshFile, "mesh file");
      Sundance::setOption("nx", nx, "number of elements in x");
      Sundance::setOption("ny", ny, "number of elements in y");
      Sundance::setOption("solver", solverFile, "name of XML file for solver");

      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      
      MeshSource mesher;
      if (meshFile != "builtin")
      {
        mesher = new ExodusMeshReader("../../../examples-tutorial/meshes/"+meshFile, meshType);
      }
      else
      {
        int npx = -1;
        int npy = -1;
        PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
        TEST_FOR_EXCEPT(npx < 1);
        TEST_FOR_EXCEPT(npy < 1);
        TEST_FOR_EXCEPT(npx * npy != np);
        mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
          0.0,  1.0, ny, npy, meshType);
      }
      Mesh mesh = mesher.getMesh();


      bool meshOK = mesh.checkConsistency(meshFile+"-check");
      if (meshOK) 
      {
        cout << "mesh is OK" << endl;
      }
      else
      {
        cout << "mesh is INCONSISTENT" << endl;
      }
      mesh.dump(meshFile+"-dump");

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left;
      CellFilter right;
      CellFilter top;
      CellFilter bottom;

      if (meshFile != "builtin")
      {
        left = edges.labeledSubset(1);
        right = edges.labeledSubset(2);
        top = edges.labeledSubset(3);
        bottom = edges.labeledSubset(4);
      }
      else
      {
        left = edges.subset(new LeftPointTest());
        right = edges.subset(new RightPointTest());
        top = edges.subset(new TopPointTest());
        bottom = edges.subset(new BottomPointTest());
      }
      
      /* Create unknown and test functions, discretized using second-order
       * Lagrange interpolants */
      BasisFamily basis = new Lagrange(1);
      Expr u = new UnknownFunction(basis, "u");
      Expr v = new TestFunction(basis, "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
      Expr one = new SundanceCore::Parameter(1.0);
      Expr exactSoln = x+2.0*y;//0.5*x*x + (1.0/3.0)*y;
      Expr eqn = Integral(interior, (grad*u)*(grad*v)  /* + one*v */, quad2);
      /* Define the Dirichlet BC */
      Expr h = new CellDiameterExpr();
      Expr bc = EssentialBC(bottom+top+left+right, v*(u-exactSoln)/h, quad4);

      Assembler::defaultVerbParams()->set<int>("global", 2);
      Assembler::defaultVerbParams()->set<int>("evaluation", 4);
      Assembler::defaultVerbParams()->set<int>("assembly loop", 2);
      
      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/" + solverFile));
#else
      ParameterXMLFileReader reader(solverFile);
#endif
      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      //      cout << "map = " << endl;
      //prob.rowMap()->print(cout);
      //cout << endl;
      Expr soln = prob.solve(solver);

      DiscreteSpace discSpace2(mesh, new Lagrange(2), vecType);
      DiscreteSpace discSpace0(mesh, new Lagrange(0), vecType);
      L2Projector proj1(discSpace2, soln-exactSoln);
      double pid = MPIComm::world().getRank();
      L2Projector proj2(discSpace0, pid);
      Expr errorDisc = proj1.project();
      Expr pidDisc = proj2.project();




      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Poisson2d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("error", new ExprFieldWrapper(errorDisc));
      w.addField("rank", new ExprFieldWrapper(pidDisc));
      w.write();

      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad4);

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
                                   derivErr*derivErr, 
                                   quad2);

      

      Expr fluxErrExpr = Integral(top, 
                                  pow(dy*(soln-exactSoln), 2),
                                  new GaussianQuadrature(2));

      

      Expr exactFluxExpr = Integral(top, 
                                    dy*exactSoln,
                                    new GaussianQuadrature(2));

      Expr numFluxExpr = Integral(top, 
                                  dy*soln,
                                  new GaussianQuadrature(2));
      //        + Integral(top+bottom, 
      //                   pow(dy*(soln-exactSoln), 2),
      //                   new GaussianQuadrature(2));


      FunctionalEvaluator errInt(mesh, errExpr);
      FunctionalEvaluator derivErrInt(mesh, derivErrExpr);

      double errorSq = errInt.evaluate();
      cout << "error norm = " << sqrt(errorSq) << endl << endl;

      double derivErrorSq = derivErrInt.evaluate();
      cout << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;

      Assembler::defaultVerbParams()->set<int>("global", 2);
      Assembler::defaultVerbParams()->set<int>("assembly loop", 2);
      double fluxErrorSq = evaluateIntegral(mesh, fluxErrExpr);
      cout << "flux error norm = " << sqrt(fluxErrorSq) << endl << endl;

      cout << "exact flux = " << evaluateIntegral(mesh, exactFluxExpr) << endl;
      cout << "numerical flux = " << evaluateIntegral(mesh, numFluxExpr) << endl;

      Sundance::passFailTest(sqrt(errorSq + derivErrorSq + fluxErrorSq), 1.0e-9);

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 

  return Sundance::testStatus();
}

#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Poisson2D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif
