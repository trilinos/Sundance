#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"


/** 
 * Solves Kepler's equation x = u(x) + e*sin(u(x))
 */

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
      const double pi = 4.0*atan(1.0);
      MeshSource mesher = new PartitionedLineMesher(0.0, pi, 20*np, meshType);
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

      /* Create coordinate function */
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

     
      /* Define the weak form */
      double e = 0.5;
      Expr eqn = Integral(interior, v*(u + e*sin(u) - x), quad);
      Expr bc;

      //      Evaluator::classVerbosity() = VerbExtreme;
      //Assembler::classVerbosity() = VerbExtreme;
      //IntegralGroup::classVerbosity() = VerbExtreme;
      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      //      F.verbosity() = VerbExtreme;
      /* Get the initial guess */
      Vector<double> x0 = F.getInitialGuess();
      
      
      /* Create an Aztec solver for solving the linear subproblems */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;
      
      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azOptions[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-13;
      
      LinearSolver<double> linSolver = new AztecSolver(azOptions,azParams);
      
 
     //  /* Set up the linear solver  */
//       ParameterList solverParams;

//       solverParams.set(LinearSolverBase<double>::verbosityParam(), 4);
//       solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
//       solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-14);

//       LinearSolver<double> linSolver = new BICGSTABSolver<double>(solverParams);


      /* Now let's create a NOX solver */

      NOX::TSF::Group grp(x0, F, linSolver);

      grp.verbosity() = VerbExtreme;

      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-15);
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


      

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  MPISession::finalize();
}
