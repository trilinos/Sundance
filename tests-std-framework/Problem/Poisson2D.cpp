#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

/** 
 * Solves the Poisson equation in 2D
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool bottomPointTest(const Point& x) {return fabs(x[1]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}
bool topPointTest(const Point& x) {return fabs(x[1]-2.0) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      VerbositySetting verb = VerbSilent;
      Assembler::classVerbosity() = verb;

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 64*np, np,
                                                         0.0, 2.0, 64, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      if (verb > VerbHigh)
        {
          FieldWriter wMesh = new VerboseFieldWriter();
          wMesh.addMesh(mesh);
          wMesh.write();
        }

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
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
      Expr one = new Parameter(1.0);
      Expr oneThird = new Parameter(1.0/3.0);
      Expr eqn = Integral(interior, (grad*v)*(grad*u)  + one*v, quad2)
        + Integral(top, -v*oneThird, quad2) 
        + Integral(right, -v*(1.5 + (1.0/3.0)*y - u), quad4);
      //        + Integral(bottom, 100.0*v*(u-0.5*x*x), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, v*(u-0.5*x*x), quad4);

      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-13;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Poisson2d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.write();

      Expr exactSoln = 0.5*x*x + (1.0/3.0)*y;

      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad4);

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
                                   derivErr*derivErr, 
                                   quad2);

      FunctionalEvaluator errInt(mesh, errExpr);
      FunctionalEvaluator derivErrInt(mesh, derivErrExpr);

      double errorSq = errInt.evaluate();
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double derivErrorSq = derivErrInt.evaluate();
      cerr << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
