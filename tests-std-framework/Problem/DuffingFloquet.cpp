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
#include "SundanceUnfoldPeriodicDF.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"


bool DuffingFloquet()
{
  int np = MPIComm::world().getNProc();
  TEUCHOS_TEST_FOR_EXCEPT(np != 1);

  const double pi = 4.0*atan(1.0);

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a periodic mesh */
  int nx = 128;

  MeshType meshType = new PeriodicMeshType1D();
  MeshSource mesher = new PeriodicLineMesher(0.0, 2.0*pi, nx, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter pts = new DimensionalCellFilter(0);
      
  CellFilter left = pts.subset(new CoordinateValueCellPredicate(0,0.0));
  CellFilter right = pts.subset(new CoordinateValueCellPredicate(0,2.0*pi));
      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u1 = new UnknownFunction(new Lagrange(1), "u1");
  Expr u2 = new UnknownFunction(new Lagrange(1), "u2");
  Expr v1 = new TestFunction(new Lagrange(1), "v1");
  Expr v2 = new TestFunction(new Lagrange(1), "v2");

  /* Create differential operator and coordinate function */
  Expr dx = new Derivative(0);
  Expr x = new CoordExpr(0);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);

  double F0 = 0.5;
  double gamma = 2.0/3.0;
  double a0 = 1.0;
  double w0 = 1.0;
  double eps = 0.5;

  Expr u1Guess = -0.75*cos(x) + 0.237*sin(x);
  Expr u2Guess = 0.237*cos(x) + 0.75*sin(x);

  DiscreteSpace discSpace(mesh, 
    List(new Lagrange(1), new Lagrange(1)),
    vecType);
  L2Projector proj(discSpace, List(u1Guess, u2Guess));
  Expr u0 = proj.project();


  Expr rhs1 = u2;
  Expr rhs2 = -w0*w0*u1 - gamma*u2 - eps*w0*w0*pow(u1,3.0)/a0/a0 
    + F0*w0*w0*sin(x);

  /* Define the weak form */
  Expr eqn = Integral(interior, 
    v1*(dx*u1 - rhs1) + v2*(dx*u2 - rhs2),
    quad);
  Expr dummyBC ; 

  NonlinearProblem prob(mesh, eqn, dummyBC, List(v1,v2), List(u1,u2), 
    u0, vecType);


  ParameterXMLFileReader reader("nox.xml");
  ParameterList solverParams = reader.getParameters();

  Out::root() << "finding periodic solution" << endl;
  NOXSolver solver(solverParams);
  prob.solve(solver);

  /* unfold the solution onto a non-periodic mesh */
      
  Expr uP = unfoldPeriodicDiscreteFunction(u0, "u_p");
  Out::root() << "uP=" << uP << endl;
      
  Mesh unfoldedMesh = DiscreteFunction::discFunc(uP)->mesh();
  DiscreteSpace unfDiscSpace = DiscreteFunction::discFunc(uP)->discreteSpace();

  FieldWriter writer = new MatlabWriter("Floquet.dat");
  writer.addMesh(unfoldedMesh);
  writer.addField("u_p[0]", new ExprFieldWrapper(uP[0]));
  writer.addField("u_p[1]", new ExprFieldWrapper(uP[1]));

  Array<Expr> a(2);
  a[0] = new Sundance::Parameter(0.0, "a1");
  a[1] = new Sundance::Parameter(0.0, "a2");


  Expr bc = EssentialBC(left, v1*(u1-uP[0]-a[0]) + v2*(u2-uP[1]-a[1]), quad);

  NonlinearProblem unfProb(unfoldedMesh, eqn, bc, 
    List(v1,v2), List(u1,u2), uP, vecType);

  unfProb.setEvalPoint(uP);

  LinearOperator<double> J = unfProb.allocateJacobian();
  Vector<double> b = J.domain().createMember();

  LinearSolver<double> linSolver
    = LinearSolverBuilder::createSolver("amesos.xml");
        
  SerialDenseMatrix<int, double> F(a.size(), a.size());

  for (int i=0; i<a.size(); i++)
  {
    Out::root() << "doing perturbed orbit #" << i << endl;
    for (int j=0; j<a.size(); j++) 
    {
      if (i==j) a[j].setParameterValue(1.0);
      else a[j].setParameterValue(0.0);
    }
        
    unfProb.computeJacobianAndFunction(J, b);
    Vector<double> w = b.copy();
    linSolver.solve(J, b, w);
    Expr w_i = new DiscreteFunction(unfDiscSpace, w);

    for (int j=0; j<a.size(); j++)
    {
      Out::root() << "postprocessing" << i << endl;

      writer.addField("w[" + Teuchos::toString(i)
        + ", " + Teuchos::toString(j) + "]", new ExprFieldWrapper(w_i[j]));
      Expr g = Integral(right, w_i[j], quad);
      F(j,i) = evaluateIntegral(unfoldedMesh, g);
    }
  }

  writer.write();

  Out::root() << "Floquet matrix = " << endl
              << F << endl;
        

  Out::root() << "doing eigenvalue analysis" << endl;
  Array<double> ew_r(a.size());
  Array<double> ew_i(a.size());
  int lWork = 6*a.size();
  Array<double> work(lWork);
  int info = 0;
  LAPACK<int, double> lapack;
  lapack.GEEV('N','N', a.size(), F.values(),
    a.size(), &(ew_r[0]), &(ew_i[0]), 0, 1, 0, 1, &(work[0]), lWork,
    &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
    std::runtime_error,
    "LAPACK GEEV returned error code =" << info);
      
  Array<double> ew(a.size());
  for (int i=0; i<a.size(); i++)
  {
    ew[i] = sqrt(ew_r[i]*ew_r[i]+ew_i[i]*ew_i[i]);
    Out::root() << setw(5) << i 
                << setw(16) << ew_r[i] 
                << setw(16) << ew_i[i] 
                << setw(16) << ew[i]
                << endl;
  }

  double err = ::fabs(ew[0] - 0.123);
  return SundanceGlobal::checkTest(err, 0.001);
}
