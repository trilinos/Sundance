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
#include "SundanceTransientStepProblem.hpp"
#include "SundanceProblemTesting.hpp"

Expr rhsF(const Expr& x, const Expr& t)
{
  return -12*pow((1 - x)*cos(t) - x*cos(t),2)
    *pow(2 + (1 - x)*x*cos(t),2) + 
    8*cos(t)*pow(2 + (1 - x)*x*cos(t),3) - (1 - x)*x*sin(t);
}



CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;});
CELL_PREDICATE(MidPointTest, {return fabs(x[0]-0.5) < 1.0e-10;});

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
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 256*np, meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);

    CellFilter right = points.subset(new RightPointTest());
    CellFilter left = points.subset(new LeftPointTest());
    CellFilter midpoint = points.subset(new MidPointTest());
      
    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    

    /* Create differential operator and coordinate function */
    Expr dx = new Derivative(0);
    Expr x = new CoordExpr(0);
    
    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);

    
    Expr tPrev = new Sundance::Parameter(0.0);
    Expr t = new Sundance::Parameter(0.0);
    Expr dt = new Sundance::Parameter(0.0);

    DiscreteSpace ds(mesh, basis, vecType);
    L2Projector proj(ds, 2.0+x*(1.0-x)*cos(t));
    Expr uPrev = proj.project();
    Expr uSoln = copyDiscreteFunction(uPrev);

    Expr eqn = Integral(interior, v*(u-uPrev) 
      + 2.0*dt*(dx*v)*(pow(u,3.0)*(dx*u) + pow(uPrev, 3.0)*(dx*uPrev))
      - 0.5*dt*v*(rhsF(x, tPrev) + rhsF(x, t)),
      quad);
    Expr bc = EssentialBC(left+right, v*(u+uPrev-4.0), quad);

    NonlinearProblem stepProb(mesh, eqn, bc, v, u, uSoln, vecType);

    TransientStepProblem prob(stepProb, tPrev, uPrev, t, uSoln, dt);

    Array<double> h;
    Array<double> err;

    for (int n=2; n<=64; n *= 2)
    {
      double t0 = 0.0;
      double deltaT = 0.5/n;
      h.append(deltaT);
      double t1 = t0 + deltaT;
      t.setParameterValue(t1);
      tPrev.setParameterValue(t0);
      dt.setParameterValue(deltaT);
      Expr uStart = proj.project();
      updateDiscreteFunction(uStart, uPrev);

      /* Use the Playa Newton-Armijo nonlinear solver */
      LinearSolver<double> linSolver = LinearSolverBuilder::createSolver("amesos.xml");
      ParameterList solverParams("NewtonArmijoSolver");
      solverParams.set("Tau Relative", 1.0e-12);
      solverParams.set("Tau Absolute", 1.0e-12);
      solverParams.set("Verbosity", 0);    
      NonlinearSolver<double> solver = new NewtonArmijoSolver<double>(solverParams, linSolver);

      for (int i=0; i<n; i++)
      {
        bool stepOK = prob.step(t0, uPrev, t1, uSoln, solver);
        t0 = t1;
        t1 += deltaT;

        FieldWriter writer = new DSVWriter("TransientNonlin-" 
          + Teuchos::toString(n) + "-" + Teuchos::toString(i) + ".dat");

        writer.addMesh(mesh);
        writer.addField("u", new ExprFieldWrapper(uSoln[0]));
        writer.write();

        updateDiscreteFunction(uSoln, uPrev);
      }
      Expr uExact = proj.project();
      Expr errExpr = Integral(interior, 
        pow(uSoln-uExact, 2),
        new GaussianQuadrature(2));
      double errorSq = evaluateIntegral(mesh, errExpr);
      Out::root() << setw(6) << n << setw(16) << t0 
                  << setw(16) << sqrt(errorSq) 
                  << setw(16) << 0.5*log(errorSq)/log(2.0) << endl;
      err.append(sqrt(errorSq));
    }
    
    double p = fitPower(h, err);

    Out::root() << "exponent = " << p << endl;

    bool pass = p > 1.95;
    Sundance::passFailTest(pass);
  }
	catch(std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}



