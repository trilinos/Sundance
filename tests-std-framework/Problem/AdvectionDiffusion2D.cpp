#include "Sundance.hpp"

/** 
 * Solves the advection-diffusion equation in 2D, with a velocity
 * field computed from a potential flow model.
 */

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}
bool bottomPointTest(const Point& x) {return fabs(x[1]) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}
bool topPointTest(const Point& x) {return fabs(x[1]-1.0) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasicSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int n = 1;
      MeshType meshType = new BasicSimplicialMeshType();
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
      int order = 2;
      Expr u = new UnknownFunction(new Lagrange(order), "u");
      Expr v = new TestFunction(new Lagrange(order), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form for the potential flow equation */
      Expr flowEqn = Integral(interior, (grad*v)*(grad*u), quad2);

      /* Define the Dirichlet BC */
      Expr flowBC = EssentialBC(bottom, v*(u-0.5*x*x), quad4)
        + EssentialBC(top, v*(u - 0.5*(x*x - 1.0)), quad4)
        + EssentialBC(left, v*(u + 0.5*y*y), quad4)
        + EssentialBC(right, v*(u - 0.5*(1.0-y*y)), quad4);

      /* We can now set up the linear problem! */
      LinearProblem flowProb(mesh, flowEqn, flowBC, v, u, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      /* solve the problem */
      Expr u0 = flowProb.solve(solver);

      /* Now set up and solve the advection-diffusion equation for r */
      Expr r = new UnknownFunction(new Lagrange(order), "u");
      Expr s = new TestFunction(new Lagrange(order), "v");

      Expr velocity = grad*u0;
      //      Expr velocity = grad*(0.5*(x*x - y*y));
      Expr adEqn = Integral(interior, (grad*s)*(grad*r), quad2)
        + Integral(interior, s*velocity*(grad*r), quad4);
        

      Expr adBC = EssentialBC(bottom, s*r, quad4)
        + EssentialBC(top, s*(r-x), quad4)
        + EssentialBC(left, s*r, quad4)
        + EssentialBC(right, s*(r-y), quad4);

      LinearProblem adProb(mesh, adEqn, adBC, s, r, vecType);
      Expr r0 = adProb.solve(solver);

      FieldWriter w = new VTKWriter("AD-2D");
      w.addMesh(mesh);
      w.addField("potential", new ExprFieldWrapper(u0[0]));
      w.addField("concentration", new ExprFieldWrapper(r0[0]));
      w.write();

      Expr exactPotential = 0.5*(x*x - y*y);
      Expr exactConcentration = x*y;

      Expr uErr = exactPotential - u0;
      Expr uErrExpr = Integral(interior, 
                              uErr*uErr,
                              quad4);

      Expr rErr = exactConcentration-r0;
      Expr rErrExpr = Integral(interior, 
                               rErr*rErr, 
                               quad4);

      FunctionalEvaluator uInt(mesh, uErrExpr);
      FunctionalEvaluator rInt(mesh, rErrExpr);

      double uErrorSq = uInt.evaluate();
      cerr << "potential error norm = " << sqrt(uErrorSq) << endl << endl;

      double rErrorSq = rInt.evaluate();
      cerr << "concentration error norm = " << sqrt(rErrorSq) << endl << endl;

      Sundance::passFailTest(uErrorSq + rErrorSq, 1.0e-11);

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
