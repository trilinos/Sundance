#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

/** 
 * Solves the Stokes equation in 2D
 */

bool leftPointTest(const Point& x) {return fabs(x[0]+1.0) < 1.0e-4;}
bool bottomPointTest(const Point& x) {return fabs(x[1]+1.0) < 1.0e-4;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-4;}
bool topPointTest(const Point& x) {return fabs(x[1]-1.0) < 1.0e-4;}

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
      int nx = 16;
      int ny = 16;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(-1.0, 1.0, nx*np, np,
                                                         -1.0, 1.0, ny, 1,
                                                         meshType);




      Mesh mesh = mesher.getMesh();
      double h = 2.0/((double) ny);

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter points = new DimensionalCellFilter(0);
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
      double beta = 0.1;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy) - p*(dx*vx+dy*vy)
                          + (h*h*beta)*(grad*q)*(grad*p) + q*(dx*ux+dy*uy),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*(1.0-y*y);
      Expr bc = EssentialBC(left, vx*ux + vy*uy, quad2)
        + EssentialBC(right, vx*ux + vy*uy, quad2)
        + EssentialBC(top, vx*(ux-y) + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


   

      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;
      //      Assembler::classVerbosity() = VerbExtreme;

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q), 
                         List(ux, uy, p), vecType);

      LinearOperator<double> A = prob.getOperator();
      ofstream ms("matrix.dat");
      A.print(ms);
      ofstream vs("vector.dat");
      Vector<double> b = prob.getRHS();
      b.print(vs);
      ofstream maps("map.dat");
      prob.rowMap()->print(maps);

      FieldWriter w1 = new VerboseFieldWriter("mesh.dat");
      w1.addMesh(mesh);
      w1.write();

      ParameterList params;
      ParameterList solverParams;
      solverParams.set("Type", "TSF");
      solverParams.set("Method", "BICGSTAB");
      solverParams.set("Max Iterations", 5000);
      solverParams.set("Restart", 100);
      solverParams.set("Tolerance", 1.0e-12);
      solverParams.set("Precond", "ILUK");
      solverParams.set("Graph Fill", 3);
      solverParams.set("Verbosity", 4);

      params.set("Linear Solver", solverParams);


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(params);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Stokes2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
