#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D
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

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int nx = 32;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx*np, np,
                                                         0.0, 1.0, nx, 1,
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
      Expr psi = new UnknownFunction(new Lagrange(1), "psi");
      Expr vPsi = new TestFunction(new Lagrange(1), "vPsi");
      Expr omega = new UnknownFunction(new Lagrange(1), "omega");
      Expr vOmega = new TestFunction(new Lagrange(1), "vOmega");

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
      Expr eqn = Integral(interior, (grad*vPsi)*(grad*psi) 
                          + (grad*vOmega)*(grad*omega) + vPsi*omega, quad2)
        + Integral(top, -1.0*vPsi, quad4);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, vOmega*psi, quad2) 
        + EssentialBC(top, vOmega*psi, quad2) 
        + EssentialBC(left, vOmega*psi, quad2) 
        + EssentialBC(right, vOmega*psi, quad2);



      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vPsi, vOmega), 
                         List(psi, omega), vecType);

      ParameterList params;
      ParameterList solverParams;
      solverParams.set("Type", "TSF");
      solverParams.set("Method", "BICGSTAB");
      solverParams.set("Max Iterations", 500);
      solverParams.set("Tolerance", 1.0e-8);
      solverParams.set("Precond", "ILUK");
      solverParams.set("Graph Fill", 1);
      solverParams.set("Verbosity", 4);

      params.set("Linear Solver", solverParams);

      XMLParameterListWriter paramWriter;
      
      cerr << "solver = " << paramWriter.toXML(params) << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(params);

    

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("VorticityStokes2d");
      w.addMesh(mesh);
      w.addField("psi", new ExprFieldWrapper(soln[0]));
      w.addField("omega", new ExprFieldWrapper(soln[1]));
      w.write();

      /* As a check, we integrate the vorticity over the domain. By 
       * Stokes' theorem this should be equal to the line integral
       * of the velocity around the boundary. */
      Expr totalVorticityExpr = Integral(interior, soln[1], quad2);
      double totalVorticity = evaluateIntegral(mesh, totalVorticityExpr);
      cerr << "total vorticity = " << totalVorticity << endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(fabs(totalVorticity-1.0), tol);
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
