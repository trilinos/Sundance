#include "Sundance.hpp"
#include "SundanceTrivialGrouper.hpp"

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
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      /* Now we set up the weak form of our equation. */
      Expr eqn = Integral(interior, v*(u*u-2.0), quad);

      /* There are no boundary conditions for this problem, so the
       * BC expression is empty */
      Expr bc;

      GrouperBase::classVerbosity() = VerbExtreme;

      /* We can now set up the nonlinear problem! */
      NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);

      /* Set up the linear solver used in solving J*delta+b = 0 */
      ParameterList solverParams;

      solverParams.set(LinearSolverBase<double>::verbosityParam(), 0);
      solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
      solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-12);

      LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);

      /* Do the nonlinear solve */
      
      Vector<double> x0 = prob.getInitialGuess();
      bool converged = false;
      for (int i=0; i<20; i++)
        {
          prob.setEvalPt(x0);
          LinearOperator<double> J = prob.getJacobian();
          Vector<double> b = prob.getFunctionValue();
          Vector<double> solnVec;
          SolverState<double> state = solver.solve(J, b, solnVec);

          cerr << "solver state = " << endl << state << endl;

          x0 = x0 - solnVec;
          cerr << "step norm = " << solnVec.norm2() << endl;

          if (solnVec.norm2() < 1.0e-14) 
            {
              cerr << "Newton's method converged!" << endl;
              converged = true;
              break;
            }
        }
      
      if (!converged) 
        {
          cerr << "FAILED TO CONVERGE!" << endl;
        }
      else
        {
          cerr << "solution is " << endl << x0 << endl;
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  MPISession::finalize();
}
