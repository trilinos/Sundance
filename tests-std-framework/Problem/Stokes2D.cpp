#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

/** 
 * Solves the Stokes equation in 2D
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
      int nx = 10;
      int ny = 10;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx*np, np,
                                                         0.0, 1.0, ny, 1,
                                                         meshType);




      Mesh mesh = mesher.getMesh();
      double h = 1.0/((double) ny);

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

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
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(1), "u_y");
      Expr vy = new TestFunction(new Lagrange(1), "v_y");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

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
      double beta = 0.025;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy) - p*(dx*vx+dy*vy)
                          - (h*h*beta)*(grad*q)*(grad*p) - q*(dx*ux+dy*uy),
                          quad2)
        + Integral(left, -(1.0-x)*vx, quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*(1.0-y*y);
      Expr bc = /* EssentialBC(left, vx*(ux-uInflow) + vy*uy , quad4)
                   + */ EssentialBC(top, vx*ux + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


   

      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;
      //      Assembler::classVerbosity() = VerbExtreme;

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q), 
                         List(ux, uy, p), vecType);


      ParameterList solverParams;

      solverParams.set(LinearSolverBase<double>::verbosityParam(), 4);
      solverParams.set(IterativeSolver<double>::maxitersParam(), 5000);
      solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-10);

      LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);

      /* Create an Aztec solver */
      // std::map<int,int> azOptions;
//       std::map<int,double> azParams;

//       azOptions[AZ_solver] = AZ_gmres;
//       azOptions[AZ_precond] = AZ_dom_decomp;
//       azOptions[AZ_subdomain_solve] = AZ_ilu;
//       azOptions[AZ_graph_fill] = 1;
//       azOptions[AZ_max_iter] = 1000;
//       azParams[AZ_tol] = 1.0e-6;

//       LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Stokes2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();

      Expr uxErr = soln[0] - uInflow;
      Expr errExpr = Integral(interior, 
                              uxErr*uxErr,
                              new GaussianQuadrature(4));

      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
