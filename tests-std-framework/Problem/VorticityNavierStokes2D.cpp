#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

/** 
 * Solves the Navier-Stokes equations on the lid-driver cavity
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool bottomPointTest(const Point& x) {return fabs(x[1]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}
bool topPointTest(const Point& x) {return fabs(x[1]-1.0) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int n=200;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n*np, np,
                                                         0.0, 1.0, n, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellPredicate rightPointFunc = new PositionalCellPredicate(rightPointTest);
      CellPredicate topPointFunc = new PositionalCellPredicate(topPointTest);
      CellPredicate bottomPointFunc = new PositionalCellPredicate(bottomPointTest);
      CellFilter left = edges.subset(leftPointFunc);
      CellFilter right = edges.subset(rightPointFunc);
      CellFilter top = edges.subset(topPointFunc);
      CellFilter bottom = edges.subset(bottomPointFunc);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr psi = new UnknownFunction(new Lagrange(1), "psi");
      Expr vPsi = new TestFunction(new Lagrange(1), "vPsi");
      Expr omega = new UnknownFunction(new Lagrange(1), "omega");
      Expr vOmega = new TestFunction(new Lagrange(1), "vOmega");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* A parameter expression for the Reynolds number */
      Expr reynolds = new Parameter(20.0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad1 = new GaussianQuadrature(1);
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      //  Assembler::classVerbosity() = VerbExtreme;

      /* Define the weak form */
      Expr psiEqn = Integral(interior, (grad*vPsi)*(grad*psi) + vPsi*omega, 
                             quad2)
        + Integral(top, 1.0*vPsi, quad1);

      Expr omegaEqn = Integral(interior, (grad*vOmega)*(grad*omega)
                               + reynolds*vOmega*((dy*psi)*dx*omega - (dx*psi)*dy*omega),
                               quad1);
      
      Expr eqn = omegaEqn + psiEqn;

      /* Define the Dirichlet BC */
      CellFilter walls = left + bottom + right;
      Expr bc = EssentialBC(walls, vOmega*psi, quad2) 
        + EssentialBC(top, vOmega*psi, quad2) ;

      BasisFamily L1 = new Lagrange(1);
      Array<BasisFamily> bases = tuple(L1, L1);
      DiscreteSpace discSpace(mesh, bases, vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      Assembler::workSetSize() = 1000;
      //      Evaluator::classVerbosity() = VerbHigh;

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, List(vPsi, vOmega),
                               List(psi, omega), u0, vecType);
      /* Get the initial guess */
      Vector<double> x0 = F.getInitialGuess();
      F.setEvalPt(x0);

      ParameterList params;
      ParameterList linSolverParams;
      linSolverParams.set("Type", "TSF");
      linSolverParams.set("Method", "BICGSTAB");
      linSolverParams.set("Max Iterations", 200);
      linSolverParams.set("Tolerance", 1.0e-14);
      linSolverParams.set("Precond", "ILUK");
      linSolverParams.set("Graph Fill", 1);
      linSolverParams.set("Verbosity", 4);

      params.set("Linear Solver", linSolverParams);


      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(params);

      NOX::TSF::Group grp(x0, F, linSolver);

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


      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("VorticityNavierStokes2d");
      w.addMesh(mesh);
      w.addField("streamfunction", new ExprFieldWrapper(u0[0]));
      w.addField("vorticity", new ExprFieldWrapper(u0[1]));
      w.write();


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
