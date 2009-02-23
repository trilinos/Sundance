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
 *
 */


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, v*(u-1.0), quad);
      Expr bc;

      /* We can now set up the linear problem! */
      cerr << "setting up linear problem" << endl;
      LinearProblem prob(mesh, eqn, bc, v, u, vecType); 


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
#else
      ParameterXMLFileReader reader("bicgstab.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      cerr << "solving problem" << endl;
      Expr soln = prob.solve(solver);

      Expr exactSoln = 1.0;

      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2),
                              new GaussianQuadrature(4));

      cerr << "setting up norm" << endl;
      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-12;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
