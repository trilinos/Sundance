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
#include "SundanceBernstein.hpp"
#include <unistd.h>
#include <sys/types.h>
/** 
 * Solves the Poisson equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      int nx = 32;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter rightPoint = points.subset(new RightPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Bernstein interpolants */
      int order = 3;
      double p = order;
      Expr u = new UnknownFunction(new Bernstein(order), "u");
      Expr v = new TestFunction(new Bernstein(order), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2*order);


      /* Define the weak form */
      Expr exactSoln = pow(1.0 + x, p);
      Expr source = p*(p-1.0)*pow(1.0+x, p-2.0);
      Expr eqn = Integral(interior, (dx*v)*(dx*u) + v*source, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*(u-exactSoln), quad)
        + EssentialBC(rightPoint, v*(u-exactSoln), quad);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType); 

      cerr << "matrix = " << endl << prob.getOperator() << endl;
      cerr << "rhs = " << endl << prob.getRHS() << endl;


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
#else
      ParameterXMLFileReader reader("bicgstab.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);


      Expr soln = prob.solve(solver);


      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2),
                              quad);

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-11;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
