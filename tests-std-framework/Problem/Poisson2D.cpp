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
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 2*np, np,
                                                         0.0, 2.0, 2, 1,
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
      Expr eqn = Integral(interior, (grad*v)*(grad*u)  + one*v, quad2)
        + Integral(top, -v/3.0, quad2) 
        + Integral(right, -v*(1.5 + (1.0/3.0)*y - u), quad4);
      //        + Integral(bottom, 100.0*v*(u-0.5*x*x), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, v*(u-0.5*x*x), quad4);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

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

      Sundance::passFailTest(errorSq + derivErrorSq, 1.0e-11);

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
