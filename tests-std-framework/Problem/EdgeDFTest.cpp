#include "Sundance.hpp"

/*
 * This program tests evaluation of more than one derivative of
 * discrete functions on submaximal cells. An error in such evaluations
 * was discovered by Paul Boggs. This test verifies the fix.
 */

bool EdgeDFTest()
{
  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

      

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedLineMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  int nx = 1;
  MeshSource mesher 
    = new PartitionedRectangleMesher(0.0, 1.0, nx, 1,
      0.0, 1.0, nx, 1, meshType);
  Mesh mesh = mesher.getMesh();

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver("amesos.xml"); 

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter edges = new DimensionalCellFilter(1);

  EvalVector::shadowOps() = true;

  Expr u = new UnknownFunction(new Lagrange(1),"u");
  Expr v = new TestFunction(new Lagrange(1), "v");

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  WatchFlag watch("watch me");
  watch.setParam("fill", 5);
  watch.setParam("integration", 5);
  watch.setParam("evaluation", 2);
  watch.setParam("discrete function evaluation", 2);
  watch.setParam("assembly loop", 5);
  watch.deactivate();

  Expr projBC;
  Expr a = x;
  Expr b = y;


  Expr projEqn1 = Integral(interior, v*(u-a), quad, watch);
  LinearProblem projProb1(mesh, projEqn1, projBC, v, u, vecType);
  Expr w0 = projProb1.solve(solver);

  Expr projEqn2 = Integral(interior, v*(u-b), quad, watch);
  LinearProblem projProb2(mesh, projEqn2, projBC, v, u, vecType);
  Expr w1 = projProb2.solve(solver);

  CellFilter cells = edges;      

  /* Create differential operator and coordinate function */
  Expr grad = gradient(2);
  Expr dx = grad[0];
  Expr tau = CellTangentExpr(2,"tau");
  Expr dt = tau*grad;

  Expr sf0 = Integral(cells, (dx*a)*(dx*a) + a, quad, watch);
  Expr sf1 = Integral(cells, (dx*w0)*(dx*w0) + w0, quad, watch);
  watch.deactivate();
  Out::os() << "============== Exact integral ================= " << endl;
  double s0 = evaluateIntegral(mesh, sf0);
  Out::os() << "============== DF integral ================= " << endl;
  double s1 = evaluateIntegral(mesh, sf1);
  Out::os() << "integral (exact) = " << s0 << endl;
  Out::os() << "integral (df)    = " << s1 << endl;
  watch.deactivate();
      
      
  double error = fabs(s1 - s0);
  double tol = 1.0e-12;
      
  return SundanceGlobal::checkTest(error, tol);
}
