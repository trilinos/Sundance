#include "Sundance.hpp"

/* 
 * test for large numbers of cell filters
 */

int main(int argc, char** argv)
{
  try
  {
    int nx = 16;
    Sundance::setOption("nx", nx, "Number of elements");

    Sundance::init(&argc, &argv);

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);
    
    Out::os() << "making vert filters" << endl;
    Array<CellFilter> X(nx+1);
    for (int i=0; i<=nx; i++)
    {
      X[i] = points.coordSubset(0, i/((double) nx));
    }
    Out::os() << "combining vert filters" << endl;
    CellFilter nodes = X[0];
    for (int i=0; i<nx; i++) 
    {
      nodes = nodes + X[i+1];
    }
    
    
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "w");
    Expr v = new TestFunction(basis, "v");

    Expr x = new CoordExpr(0);

    QuadratureFamily quad2 = new GaussianQuadrature(2);

    Expr eqn 
      = Integral(interior, v*(u-x), quad2) + Integral(nodes, v*(u-x), quad2);

    Expr bc;

    Out::os() << "forming prob" << endl;

    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    LinearSolver<double> linSolver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    Out::os() << "solving prob" << endl;
    Expr soln = prob.solve(linSolver);
    Out::os() << "done solving prob" << endl;
  }
	catch(exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return 0;
}

