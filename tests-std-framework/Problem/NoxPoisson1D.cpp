#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "NOX_TSF_StatusTestBuilder.H"

#include "TSFNOXSolver.H"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"

/** 
 * Solves the Poisson equation in 2D using the NOX solver. 
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}

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
      CellFilter points = new DimensionalCellFilter(0);
      CellPredicate leftPointFunc = new PositionalCellPredicate(leftPointTest);
      CellFilter leftPoint = points.subset(leftPointFunc);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, -(dx*v)*(dx*u) - 2.0*v, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*u, quad);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");


      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      
      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/nox.xml");
      ParameterList noxParams = reader.getParameters();

      cerr << "solver params = " << noxParams << endl;

      NOXSolver solver(noxParams, F);

      solver.solve();

      Expr exactSoln = x*(x-2.0);

      Expr errExpr = Integral(interior, 
                              pow(u0-exactSoln, 2),
                              new GaussianQuadrature(4));

      Expr derivErrExpr = Integral(interior, 
                                   pow(dx*(u0-exactSoln), 2),
                                   new GaussianQuadrature(2));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
      cerr << "deriv error norm = " << sqrt(derivErrorSq) << endl << endl;
      
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  MPISession::finalize();
}

