/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"

/* ------------------------------------------------------------------------ 
 *
 * This program solves a nonlinear initial-boundary-value problem in 
 * one spatial dimension using Galerkin finite elements in space and
 * Crank-Nicolson stepping in time. 
 *
 * ------------------------------------------------------------------------ */

const double pi = 4.0*atan(1.0);

Expr force(const double& epsilon, const Expr& x, const Expr& t)
{
  Expr rtn = -8.0*epsilon*cos(2.0*pi*t)*
        pow(1.0 + pow(x,2.0)*epsilon*cos(2.0*pi*t),2.0)*
        (1.0 + 7.0*pow(x,2.0)*epsilon*cos(2.0*pi*t)) - 
        2.0*pi*pow(x,2.0)*epsilon*sin(2.0*pi*t);

  return rtn;
}




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
      int nx = 32;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter rightPoint = points.subset(new RightPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      BasisFamily bas = new Lagrange(1);
      Expr u = new UnknownFunction(bas, "u");
      Expr v = new TestFunction(bas, "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      double epsilon = 0.25;
      /* The initial profile is u(x,0)=1 + epsilon x^2. 
       * Project this onto a discrete function */
      DiscreteSpace discSpace(mesh, bas, vecType);
      L2Projector projector(discSpace, 1.0 + epsilon*x*x);
      Expr uPrev = projector.project();

      Expr uStart = copyDiscreteFunction(uPrev, "uStart");
      Expr uNewt = copyDiscreteFunction(uPrev, "uNewt");

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      int nSteps = 32;
      double dt = 1.0/((double) nSteps);
      Expr t = new Sundance::Parameter(0.0);
      Expr tPrev = new Sundance::Parameter(0.0);

      /* Define the weak form, semidiscretized in time */
      Expr eqn = Integral(interior, v*(u-uPrev) 
        + dt/2.0*(dx*v)*((dx*pow(u, 4.0))+(dx*pow(uPrev, 4.0)))
        - dt/2.0*v*(force(epsilon, x, t)+force(epsilon, x, tPrev)), quad); 
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*(u-1.0), quad) 
        + EssentialBC(rightPoint, v*(u - 1.0 - epsilon*cos(2.0*pi*t)),quad);

      /* We can now set up the linear problem! */
      NonlinearProblem prob(mesh, eqn, bc, v, u, uNewt, vecType); 

      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver("amesos.xml");


      LinearOperator<double> J = prob.allocateJacobian();
      Vector<double> resid = J.range().createMember();
      Vector<double> newtonStep = J.domain().createMember();

      int maxNewtIters = 10;
      double newtTol = 1.0e-12;

      FieldWriter w = new MatlabWriter("transientNonlinear1D-0.dat"); 
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(uPrev[0]));
      w.addField("uNewt", new ExprFieldWrapper(uNewt[0]));
      w.write();

      /* loop over timesteps */
      for (int i=0; i<nSteps; i++)
      {
        Out::root() << "timestep #" << i << endl;
        t.setParameterValue((i+1)*dt);
        tPrev.setParameterValue(i*dt);

        /* loop over Newton steps */
        bool newtonConverged = false;

        for (int j=0; j<maxNewtIters; j++)
        {
          prob.setInitialGuess(uNewt);
          prob.computeJacobianAndFunction(J, resid);
          SolverState<double> solveState 
            = linSolver.solve(J, resid, newtonStep);
          TEST_FOR_EXCEPTION(solveState.finalState() != SolveConverged,
            std::runtime_error,
            "linear solve failed!");
          addVecToDiscreteFunction(uNewt, -1.0*newtonStep);
          double newtStepNorm = newtonStep.norm2();
          Out::root() << "|newt step| = " << newtStepNorm << endl;
          if (newtStepNorm < newtTol) 
          {
            newtonConverged = true;
            break;
          }
        }

        TEST_FOR_EXCEPTION(!newtonConverged, std::runtime_error,
          "Newton's method failed to converged after " << maxNewtIters
          << " iterations");

        updateDiscreteFunction(uNewt, uPrev);

        FieldWriter writer = new MatlabWriter("transientNonlinear1D-" 
          + Teuchos::toString(i+1) + ".dat");
        writer.addMesh(mesh);
        writer.addField("u", new ExprFieldWrapper(uPrev[0]));
        writer.write();
      }

      double err = L2Norm(mesh, interior, uPrev - uStart, quad);
      Out::root() << "error = " << err << endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(err, tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
