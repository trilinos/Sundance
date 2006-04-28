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


/** 
 * Solves the radiation diffusion equation in 1D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, void** argv)
{
 
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 20, 2, 
                                                         0.0, 1.0, 20, 2, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter rightPoint = points.subset(new RightPointTest());
      CellFilter leftPoint = points.subset(new LeftPointTest());
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);

      /* Create a discrete space, and discretize the function 1.0+x on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector projector(discSpace, 1.0+x);
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(8);

     
      /* Define the weak form */
      Expr eqn = Integral(interior, u*u*u*(grad*v)*(grad*u), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*(u-(x+1.0)), quad)
        + EssentialBC(rightPoint, v*(u-(x+1.0)), quad); 

  

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      
      ParameterXMLFileReader reader("../../../examples-tutorial/data/nox-aztec.xml");
      ParameterList noxParams = reader.getParameters();

      cerr << "solver params = " << noxParams << endl;

      NOXSolver solver(noxParams, F);

      solver.solve();

      /* check solution */
      Expr exactSoln = pow(15.0*x + 1.0, 0.25);
      
      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(8));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      

      

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
