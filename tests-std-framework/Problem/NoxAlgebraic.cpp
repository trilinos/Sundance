/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"

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

      /* We can now set up the nonlinear problem! */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);

      /* Get the initial guess */
      Vector<double> x0 = F.getInitialGuess();


      /* Create an Aztec solver for solving the linear subproblems */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-13;

      LinearSolver<double> linSolver = new AztecSolver(azOptions,azParams);



      /* Now let's create a NOX solver */

      NOX::TSF::Group grp(x0, F, linSolver);

      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-10);
      NOX::StatusTest::MaxIters statusTestB(20);
      NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

      // Create the list of solver parameters
      NOX::Parameter::List solverParameters;

      // Set the solver (this is the default)
      solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

      // Create the line search parameters sublist
      NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

      // Set the line search method
      lineSearchParameters.setParameter("Method","More'-Thuente");

      // Create the solver
      NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = solver.solve();

      // Print the answer
      cout << "\n" << "-- Parameter List From Solver --" << "\n";
      solver.getParameterList().print(cout);

      // Get the answer
      grp = solver.getSolutionGroup();

      // Print the answer
      cout << "\n" << "-- Final Solution From Solver --" << "\n";
      grp.print();


      cerr << "|F| is " << endl << grp.getNormF() << endl;
      double tol = 1.0e-10;
      Sundance::passFailTest(grp.getNormF(), tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  MPISession::finalize();
}
