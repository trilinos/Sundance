#include "Sundance.hpp"


/** 
 * Solves the radiation diffusion equation in 1D
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}

int main(int argc, void** argv)
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
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 100*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellPredicate rightPointFunc = new PositionalCellPredicate(rightPointTest);
      CellFilter leftPoint = points.subset(leftPointFunc);
      CellFilter rightPoint = points.subset(rightPointFunc);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* Create a discrete space, and discretize the function 1.0+x on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector projector(discSpace, 1.0+x);
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(8);

     
      /* Define the weak form */
      Expr eqn = Integral(interior, u*u*u*(dx*v)*(dx*u), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*(u-1.0), quad)
        + EssentialBC(rightPoint, v*(u-2.0), quad); 

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/nox.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams, F);

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = solver.solve();

      /* check solution */
      Expr exactSoln = pow(15.0*x + 1.0, 0.25);
      
      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(8));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      

      double tol = 1.0e-4;
      Sundance::passFailTest(errorSq, tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
