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
#include "SundanceFunctional.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"

using SundanceCore::List;
/** 
 *
 */

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
      int nx = 8;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx*np, np,
                                                         0.0, 1.0, nx, 1,
                                                         meshType);

      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      
      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Compute an integral over a fixed integrand */
      const double pi = 4.0*atan(1.0);

//       Expr I1 = Integral(interior, x*sin(pi*x), quad4);
//       double f1 = evaluateIntegral(mesh, I1);
//       cerr << "integral of x sin(pi*x) = " << f1 << endl;
//       double I1Exact = 1.0/pi;
//       cerr << "exact: " << I1Exact << endl;

//       double error = fabs(f1 - I1Exact);
//       cerr << "error = " << fabs(f1 - I1Exact) << endl;

//       Expr I2 = Integral(interior, x*x*sin(pi*x), quad4);
//       double f2 = evaluateIntegral(mesh, I2);
//       cerr << "integral of x^2 sin(pi*x) = " << f2 << endl;
//       double I2Exact = (1.0 - 4.0/pi/pi)/pi;
//       cerr << "exact: " << I2Exact << endl;

//       error = max(error, fabs(f2 - I2Exact));
//       cerr << "error = " << fabs(f2 - I2Exact) << endl;

//       Expr I3 = Integral(interior, sin(pi*x)*sin(pi*x), quad4);
//       double f3 = evaluateIntegral(mesh, I3);
//       cerr << "integral of sin^2(pi*x) = " << f3 << endl;
//       double I3Exact = 0.5;
//       cerr << "exact: " << I3Exact << endl;

//       error = max(error, fabs(f3 - I3Exact));
//       cerr << "error = " << fabs(f3 - I3Exact) << endl;

//       Expr I4 = Integral(interior, x*x*(pi-x)*(pi-x), quad4);
//       double f4 = evaluateIntegral(mesh, I4);
//       cerr << "integral of x^2 (pi-x)^2 = " << f4 << endl;
//       double I4Exact = pi*pi/3.0 - pi/2 + 1.0/5.0;
//       cerr << "exact: " << I4Exact << endl;

//       error = max(error, fabs(f4 - I4Exact));
//       cerr << "error = " << fabs(f4 - I4Exact) << endl;


      /* now compute a functional at a particular value of a field */
      Expr alpha = new UnknownFunction(new Lagrange(2), "u");

      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      L2Projector projector(discSpace, x*(pi-x));
      Expr alpha0 = projector.project();

      Expr g = Integral(interior, 0.5*pow(alpha-sin(pi*x), 2.0) , quad4);

      //   double gExact = 0.5*(-2.0*pi*I1Exact + I3Exact + 2.0*I2Exact + I4Exact);
      
      Functional G(mesh, g, vecType);

      FunctionalEvaluator gEval = G.evaluator(alpha, alpha0);
      double gVal = gEval.evaluate();
      cerr << "integral value = " << gVal << endl;
      //      cerr << "exact value = " << gExact << endl;
      //      error = max(error, fabs(gVal - gExact));
      //      cerr << "error = " << fabs(gVal - gExact) << endl;

      /* now compute the derivative of a functional wrt a field variable */

      Expr dGdAlpha = gEval.evalGradient(gVal);
      cerr << "integral value = " << gVal << endl;
      //      error = max(error, fabs(gVal - gExact));
      //      cerr << "error = " << fabs(gVal - gExact) << endl;

      double h = 1.0e-2;
      double diffErr = gEval.fdGradientCheck(h);

  //     error = max(error, fabs(diffErr));
//       cerr << "error = " << error << endl;

      double tol = 1.0e-8;
      Sundance::passFailTest(diffErr, tol);


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
