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

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;});
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;});
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-2.0) < 1.0e-10;});


int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int nx = 40;
      int ny = 40;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
                                                         0.0, 2.0, ny, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());

      
      /* Create unknown and test functions, discretized using second-order
       * Lagrange interpolants */
      BasisFamily basis = new Lagrange(2);
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
      Expr one = new Parameter(1.0);
      Expr eqn = Integral(interior, (dx*u)*(dx*v) + (dy*u)*(dy*v)  + one*v, quad2)
        + Integral(top, -v/3.0, quad2) 
        + Integral(right, -v*(1.5 + (1.0/3.0)*y - u), quad4);
      //        + Integral(bottom, 100.0*v*(u-0.5*x*x), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, v*(u-0.5*x*x), quad4);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/aztec.xml");
      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      //      cout << "map = " << endl;
      //prob.rowMap()->print(cout);
      //cout << endl;
      Expr soln = prob.solve(solver);

      Expr exactSoln = 0.5*x*x + (1.0/3.0)*y;

      /* Write the field in VTK format */
      // FieldWriter w = new VTKWriter("Poisson2d");
//       w.addMesh(mesh);
//       w.addField("soln", new ExprFieldWrapper(soln[0]));
//       w.write();

      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad4);

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
                                   derivErr*derivErr, 
                                   quad2);

      FunctionalEvaluator errInt(mesh, errExpr);
      FunctionalEvaluator derivErrInt(mesh, derivErrExpr);

      double errorSq = errInt.evaluate();
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double derivErrorSq = derivErrInt.evaluate();
      cerr << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;

      Sundance::passFailTest(errorSq + derivErrorSq, 1.0e-11);

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
