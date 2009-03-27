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

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int nx = 4;
      int ny = 4;
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEST_FOR_EXCEPT(npx < 1);
      TEST_FOR_EXCEPT(npy < 1);
      TEST_FOR_EXCEPT(npx * npy != np);
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                                         0.0, 1.0, ny, npy,
                                                         meshType);

      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      
      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad4 = new GaussianQuadrature(6);

      /* Compute an integral over a fixed integrand */
      const double pi = 4.0*atan(1.0);

      Expr I0 = Integral(interior, x, quad4);
      double f0 = evaluateIntegral(mesh, I0);
      cout << "integral of x = " << f0 << endl;
      double I0Exact = 0.5;
      cout << "exact: " << I0Exact << endl;

      double error = fabs(f0 - I0Exact);
      cout << "error = " << fabs(f0 - I0Exact) << endl;

      Expr I1 = Integral(interior, x*sin(pi*x), quad4);
      double f1 = evaluateIntegral(mesh, I1);
      cout << "integral of x sin(pi*x) = " << f1 << endl;
      double I1Exact = 1.0/pi;
      cout << "exact: " << I1Exact << endl;

      error = max(error, fabs(f1 - I1Exact));
      cout << "error = " << fabs(f1 - I1Exact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I2 = Integral(interior, x*x*sin(pi*x), quad4);
      double f2 = evaluateIntegral(mesh, I2);
      cout << "integral of x^2 sin(pi*x) = " << f2 << endl;
      double I2Exact = (1.0 - 4.0/pi/pi)/pi;
      cout << "exact: " << I2Exact << endl;

      error = max(error, fabs(f2 - I2Exact));
      cout << "error = " << fabs(f2 - I2Exact) << endl; 

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I3 = Integral(interior, sin(pi*x)*sin(pi*x), quad4);
      double f3 = evaluateIntegral(mesh, I3);
      cout << "integral of sin^2(pi*x) = " << f3 << endl;
      double I3Exact = 0.5;
      cout << "exact: " << I3Exact << endl;

      error = max(error, fabs(f3 - I3Exact));
      cout << "error = " << fabs(f3 - I3Exact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I4 = Integral(interior, x*x*(pi-x)*(pi-x), quad4);
      double f4 = evaluateIntegral(mesh, I4);
      cout << "integral of x^2 (pi-x)^2 = " << f4 << endl;
      double I4Exact = pi*pi/3.0 - pi/2 + 1.0/5.0;
      cout << "exact: " << I4Exact << endl;

      error = max(error, fabs(f4 - I4Exact));
      cout << "error = " << fabs(f4 - I4Exact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();


      Expr I5 = Integral(interior, 1.0, quad4);
      double f5 = evaluateIntegral(mesh, I5);
      cout << "integral of 1.0 = " << f5 << endl;
      double I5Exact = 1.0;
      cout << "exact: " << I5Exact << endl;

      error = max(error, fabs(f5 - I5Exact));
      cout << "error = " << fabs(f5 - I5Exact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();


      /* now compute a functional at a particular value of a field */
      Expr alpha = new UnknownFunction(new Lagrange(2), "alpha");
      Expr beta = new TestFunction(new Lagrange(2), "beta");

      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      L2Projector projector(discSpace, x*(pi-x));
      Expr alpha0 = projector.project();

      cout << "computing Integral(0.5*pow(alpha-sin(pi*x), 2)) at alpha0=x*(pi-x)"
           << endl;
      //#define BLAHBLAH 1
#ifdef LOUD
      verbosity<Evaluator>() = VerbExtreme;
      verbosity<SparsitySuperset>() = VerbExtreme;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<Assembler>() = VerbExtreme;
      verbosity<ElementIntegral>() = VerbExtreme;
      EvalVector::shadowOps() = true;
#endif   
      Expr g = Integral(interior, 0.5*pow(alpha-sin(pi*x), 2.0) , quad4);
      //    Expr g = Integral(interior, 0.5*(alpha-sin(pi*x))*(alpha-sin(pi*x)) , quad4);
      //Expr g = Integral(interior, sin(alpha), quad4);
      //Expr dg = Integral(interior, alpha*beta+cos(alpha0)*beta, quad4);
      Expr dg = Integral(interior, alpha*beta 
                         + (alpha0-sin(pi*x))*beta , quad4);
      Expr bc;
      LinearProblem phony(mesh, dg, bc, beta, alpha, vecType);
      Vector<double> dgVec = phony.getRHS();

      double gExact = 0.5*(-2.0*pi*I1Exact + I3Exact + 2.0*I2Exact + I4Exact);

      //#define BLAHBLAH 1
#ifdef BLAHBLAH
      verbosity<Evaluator>() = VerbExtreme;
      verbosity<SparsitySuperset>() = VerbExtreme;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<Assembler>() = VerbExtreme;
      verbosity<ElementIntegral>() = VerbExtreme;
      EvalVector::shadowOps() = true;
#endif   

      Functional G(mesh, g, vecType);

      FunctionalEvaluator gEval = G.evaluator(alpha, alpha0);

      cout << "computing function value: " << endl;
      double gVal = gEval.evaluate();
      cout << "integral value = " << gVal << endl;
      cout << "exact value = " << gExact << endl;
      error = max(error, fabs(gVal - gExact));
      cout << "error = " << fabs(gVal - gExact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();
      
      /* now compute the derivative of a functional wrt a field variable */
 
      cout << "computing function value and gradient together: " << endl;
      Expr dGdAlpha = gEval.evalGradient(gVal);
      cout << "getting vector " << endl;
      Vector<double> dgNumVec 
        = DiscreteFunction::discFunc(dGdAlpha)->getVector();
      Vector<double> dgDiff = dgVec - dgNumVec;
      cout << "grad diff = " << endl << dgDiff.norm2() << endl; 

      cout << "integral value = " << gVal << endl;
      error = max(error, fabs(gVal - gExact));
      cout << "error = " << fabs(gVal - gExact) << endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();
      cout << "*********************** FD check ***************************** " << endl;
      double h = 1.0e-2;
      double diffErr = gEval.fdGradientCheck(h);

      error = max(error, fabs(diffErr));
      cout << "max error = " << error << endl;

      double tol = 1.0e-8;
      Sundance::passFailTest(error, tol);


    }
	catch(exception& e)
		{
      cout << e.what() << endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
