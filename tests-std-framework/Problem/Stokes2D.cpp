/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
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
      int nx = 64;
      int ny = 64;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(-1.0, 1.0, nx*np, np,
                                                         -1.0, 1.0, ny, 1,
                                                         meshType);




      Mesh mesh = mesher.getMesh();
      double h0 = 2.0/((double) ny);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr h = new CellDiameterExpr();

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

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      double beta = 0.02;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + (h*h*beta)*(grad*q)*(grad*p) + q*(dx*ux+dy*uy),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*(1.0-y*y);
      Expr bc = EssentialBC(left, vx*(ux-uInflow) + vy*uy, quad2)
        + EssentialBC(top, vx*ux + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


   

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q), 
                         List(ux, uy, p), vecType);

     

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);


      Expr exactUx = uInflow;
      Expr exactUy = 0.0;
      Expr errX = exactUx - soln[0];
      Expr errY = exactUy - soln[1];

      Expr errXExpr = Integral(interior, 
                              errX*errX,
                              new GaussianQuadrature(6));

      Expr errYExpr = Integral(interior, 
                              errY*errY,
                              new GaussianQuadrature(6));

      FunctionalEvaluator errXInt(mesh, errXExpr);
      FunctionalEvaluator errYInt(mesh, errYExpr);

      double errorXSq = errXInt.evaluate();
      double errorYSq = errYInt.evaluate();
      cerr << "error norm |u_x - u_x(0)| = " << sqrt(errorXSq) << endl << endl;
      cerr << "error norm |u_y - u_y(0)| = " << sqrt(errorYSq) << endl << endl;

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Stokes2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();

      

      double tol = 1.0e-4;
      Sundance::passFailTest(sqrt(errorXSq+errorYSq), tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
