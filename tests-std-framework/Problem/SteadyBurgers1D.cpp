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
#include "SundanceUnknownParameter.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionData.hpp"


/** 
 * Solves the steady Burgers equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 400*np, meshType);
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
      Expr x = new CoordExpr(0);

      /* Parameters */
      Expr a = new UnknownParameter("a");
      Expr b = new UnknownParameter("b");
      Expr c = new UnknownParameter("c");
      Expr a0 = new SundanceCore::Parameter(-2.0);
      Expr b0 = new SundanceCore::Parameter(1.0);
      Expr c0 = new SundanceCore::Parameter(10.0);

      Expr uLeft = -2.0*c*b/a;
      Expr uRight = -2.0*c*b/(a+b);
      Expr exactSoln = -2.0*c0/(a0 + b0*x);

      Expr den = a0+b0*x;
      Expr sa = 2.0*b0*c0/den/den;
      Expr sb = -2.0*c0*a0/den/den;
      Expr sc = -2.0*b0/den;

      /* Create a discrete space, and discretize the function 1.0+x on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector projector(discSpace, x);
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

     
      /* Define the weak form */
      Expr eqn = Integral(interior, c*(dx*u)*(dx*v) + v*u*(dx*u), quad);
      /* Define the Dirichlet BC */

      WatchFlag watch("watch eqn");

      Expr bc = EssentialBC(leftPoint, v*(u-uLeft), quad, watch)
        + EssentialBC(rightPoint, v*(u-uRight), quad); 

      /* Create a TSF NonlinearOperator object */
      NonlinearProblem prob(mesh, eqn, bc, v, u, u0, 
        List(a,b,c), List(a0,b0,c0), vecType);

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);
      LinearSolver<double> linSolver = solver.linSolver();

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = prob.solve(solver);
      TEST_FOR_EXCEPTION(status != NOX::StatusTest::Converged,
        runtime_error, "solve failed");


      /* compute senstivities */
      Expr sens = prob.computeSensitivities(linSolver);

      /* Write the field in ASCII format */
      FieldWriter w = new MatlabWriter("Burgers1DSoln");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0));
      for (int i=0; i<sens.size(); i++)
      {
        w.addField("sens_" + Teuchos::toString(i), 
          new ExprFieldWrapper(sens[i]));
      }
      w.write();


      /* check solution */
      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(8));
      Expr errExprA = Integral(interior, 
                              pow(sens[0]-sa, 2),
                              new GaussianQuadrature(8));
      Expr errExprB = Integral(interior, 
                              pow(sens[1]-sb, 2),
                              new GaussianQuadrature(8));
      Expr errExprC = Integral(interior, 
                              pow(sens[2]-sc, 2),
                              new GaussianQuadrature(8));

      double errorSq0 = evaluateIntegral(mesh, errExpr);
      cerr << "soln error norm = " << sqrt(errorSq0) << endl << endl;

      double errorSqA = evaluateIntegral(mesh, errExprA);
      cerr << "sens A error norm = " << sqrt(errorSqA) << endl << endl;

      double errorSqB = evaluateIntegral(mesh, errExprB);
      cerr << "sens B error norm = " << sqrt(errorSqB) << endl << endl;

      double errorSqC = evaluateIntegral(mesh, errExprC);
      cerr << "sens C error norm = " << sqrt(errorSqC) << endl << endl;


      double error = sqrt(errorSq0 + errorSqA + errorSqB + errorSqC);
      
      double tol = 1.0e-4;
      Sundance::passFailTest(error, tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
