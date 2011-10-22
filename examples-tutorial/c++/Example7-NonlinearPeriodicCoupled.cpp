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
#include "SundancePeriodicLineMesher.hpp"
#include "SundancePeriodicMeshType1D.hpp"

/** 
 * Solves the equation
 * 
 * \f[ u'' + 2 u' + u^2 = \sin^2(2x+1) - 4 \sin(2x+1) + 4 \cos(2x+1) \f]
 *
 * with periodic BC on the interval \f$ (0, 2\pi) \f$.
 */

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      TEUCHOS_TEST_FOR_EXCEPT(np != 1);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a periodic mesh */
      int nx = 32;
      const double pi = 4.0*atan(1.0);
      MeshType meshType = new PeriodicMeshType1D();
      MeshSource mesher = new PeriodicLineMesher(0.0, 2.0*pi, nx, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u1 = new UnknownFunction(new Lagrange(1), "u1");
      Expr u2 = new UnknownFunction(new Lagrange(1), "u2");
      Expr v1 = new TestFunction(new Lagrange(1), "v1");
      Expr v2 = new TestFunction(new Lagrange(1), "v2");

      Expr t1 = new TestFunction(new Lagrange(1), "t1");
      Expr t2 = new TestFunction(new Lagrange(1), "t2");
      Expr s1 = new UnknownFunction(new Lagrange(1), "s1");
      Expr s2 = new UnknownFunction(new Lagrange(1), "s2");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      double b1 = 1.0;
      double b2 = 1.0;
      double phi1 = 0.0;
      double phi2 = pi/4.0;

      Expr P1 = b1 + cos(x-phi1);
      Expr P2 = b2 + cos(x-phi2);
      Expr Q1 = exp(-sin(x-phi1));
      Expr Q2 = exp(-sin(x-phi2));

      Expr r1 = (2.0*u2 - u1)/3.0;
      Expr r2 = (2.0*u1 - u2)/3.0;
      Expr H1 = r1*(P1 - Q1*r1);
      Expr H2 = r2*(P2 - Q2*r2);

      Expr rhs1 = (H1+2.0*H2);
      Expr rhs2 = (2.0*H1+H2);
      Out::root() << "rhs1 = " << rhs1 << endl;
      Out::root() << "rhs2 = " << rhs2 << endl;

      /* Define the weak form */
      Expr eqn = Integral(interior, 
        v1*(dx*u1 - rhs1) + v2*(dx*u2 - rhs2)
        + t1*(dx*s1 - s1*(P1-s1*Q1)) + t2*(dx*s2 - s2*(P2-s2*Q2)),
        quad);
      Expr bc ; // no explicit BC needed

      /* We can now set up the linear problem! */

      DiscreteSpace discSpace(mesh, 
        List(new Lagrange(1), new Lagrange(1),
          new Lagrange(1), new Lagrange(1)),
        vecType);
      Expr u0 = new DiscreteFunction(discSpace, 2.0);

      NonlinearProblem prob(mesh, eqn, bc, List(v1,v2,t1,t2), List(u1,u2,s1,s2), 
        u0, vecType);


      ParameterXMLFileReader reader("nox.xml");
      ParameterList solverParams = reader.getParameters();

      NOXSolver solver(solverParams);
      prob.solve(solver);

      Expr uEx1 = b1/Q1 + 2.0*b2/Q2;
      Expr uEx2 = 2.0*b1/Q1 + b2/Q2;



      L2Projector proj(discSpace, List(uEx1, uEx2, b1/Q1, b2/Q2));
      Expr uEx0 = proj.project();

      Vector<double> v0 = getDiscreteFunctionVector(u0);
      Vector<double> vEx0 = getDiscreteFunctionVector(uEx0);
      Out::root() << setw(5) << "i" << setw(16) << "x"
                  << setw(16) << "u1" << setw(16) << "u1Exact"
                  << setw(16) << "u2" << setw(16) << "u2Exact"
                  << setw(16) << "s1" << setw(16) << "s1Exact" 
                  << setw(16) << "s2" << setw(16) << "s2Exact" << endl;
      for (int i=0; i<nx; i++)
      {
        if (i>0 && i%5==0) Out::root() << endl;
        Out::root() <<setw(5) << i <<  setw(16) << i*(2.0*pi/nx) 
                    << setw(16) << v0[4*i] << setw(16) << vEx0[4*i] 
                    << setw(16) << v0[4*i+1] << setw(16) << vEx0[4*i+1] 
                    << setw(16) << v0[4*i+2] << setw(16) << vEx0[4*i+2] 
                    << setw(16) << v0[4*i+3] << setw(16) << vEx0[4*i+3] 
                    << endl;
      }
      
      Expr uErrExpr = Integral(interior, 
        pow(uEx0[0]-u0[0],2), new GaussianQuadrature(6));
      
      FunctionalEvaluator uErrInt(mesh, uErrExpr);

      double uErrorSq = uErrInt.evaluate();
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      double tol = 1.0e-3;
      Sundance::passFailTest(sqrt(uErrorSq), tol);

    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
