/** @HEADER@ */
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

using Sundance::List;

using namespace Sundance;
using namespace Playa;


/**
 *
 * This example shows how to sweep over a changing parameter without creating a new
 * expression and problem for each parameter value. 
 *
 * Similar logic can be used for continuation methods or parameter space sampling.
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
      int nx = 32;

      MeshSource mesher = new PartitionedRectangleMesher(
        0.0, 1.0, nx,
        0.0, 1.0, nx,
        meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      /* Make cell filters for the east and west boundaries */
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter west = edges.coordSubset(0, 0.0);
      CellFilter east = edges.coordSubset(0, 1.0);

      /* Create unknown and test functions */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr x = new CoordExpr(0);
      Expr grad = gradient(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      /* Define the parameter */
      Expr xi = new Sundance::Parameter(0.0);

      /* Construct a forcing term to provide an exact solution for
       * validation purposes. This will involve the parameter. */
      Expr uEx = x*(1.0-x)*(1.0+xi*exp(x));
      Expr f = -(-20 - exp(x)*xi*(1 + 32*x + 10*x*x + 
          exp(x)*(-1 + 2*x*(2 + x))*xi))/10.0;

      /* Define the weak form, using the parameter expression. This weak form
       * can be used for all parameter values. */
      Expr eqn = Integral(interior, 
        (1.0 + 0.1*xi*exp(x))*(grad*v)*(grad*u) - f*v, quad);

      /* Define the Dirichlet BC */
      Expr h = new CellDiameterExpr();
      Expr bc = EssentialBC(east + west, v*u/h, quad);

      /* We can now set up the linear problem. This can be reused
       * for different parameter values. */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      /* make a projector for the exact solution. Just like the
       * problem, this can be reused for different parameter values. */
      DiscreteSpace ds(mesh, new Lagrange(1), vecType);
      L2Projector proj(ds, uEx);

      /* Get the solver and declare variables for the results */
      LinearSolver<double> solver = LinearSolverBuilder::createSolver("aztec-ml.xml");
      Expr soln;
      SolverState<double> state;

      /* Set up the sweep from xi=0 to xi=xiMax in nSteps steps. */
      int nSteps = 10;
      double xiMax = 2.0;
      
      /* Make an array in which to keep the observed errors */
      Array<double> err(nSteps);

      /* Do the sweep */
      for (int n=0; n<nSteps; n++)
      {
        /* Update the parameter value */
        double xiVal = xiMax*n/(nSteps - 1.0);
        xi.setParameterValue(xiVal);
        Out::root() << "step n=" << n << " of " << nSteps << " xi=" << xiVal;

        /* Solve the problem. The updated parameter value is automatically used. */
        state = prob.solve(solver, soln);

        TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
          std::runtime_error,
          "solve failed!");

        /* Project the exact solution onto a discrrete space for viz. The updated
         * parameter value is automatically used. */
        Expr uEx0 = proj.project();

        /* Write the approximate and exact solutions for viz */
        FieldWriter w = new VTKWriter("ParameterSweep-" + Teuchos::toString(n));
        w.addMesh(mesh);
        w.addField("u", new ExprFieldWrapper(soln[0]));
        w.addField("uEx", new ExprFieldWrapper(uEx0[0]));
        w.write();

        /* Compute the L2 norm of the error */
        err[n] = L2Norm(mesh, interior, soln-uEx, quad);
        Out::root() << " L2 error = " << err[n] << endl;
      } 

      /* The errors are O(h^2), so use that to set a tolerance */
      double hVal = 1.0/(nx-1.0);
      double fudge = 2.0;
      double tol = fudge*hVal*hVal;

      /* Find the max error over all parameter values */
      double maxErr = *std::max_element(err.begin(), err.end());

      /* Check the error */
      Sundance::passFailTest(maxErr, tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); 
  return Sundance::testStatus(); 
}
