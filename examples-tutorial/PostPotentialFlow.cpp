#include "Sundance.hpp"

/** 
 * Solves the Laplace equation for potential flow past a dome in 
 * a wind tunnel.
 */


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
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/finePost.ncdf", meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter boundary = new BoundaryCellFilter();
      CellFilter in = boundary.labeledSubset(3);
      CellFilter out = boundary.labeledSubset(5);
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr phi = new UnknownFunction(new Lagrange(1), "u");
      Expr phiHat = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate functions */
      Expr x = new CoordExpr(0);
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      
      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      double L = 1.0;
      /* Define the weak form */
      Expr eqn = Integral(interior, (grad*phiHat)*(grad*phi), quad2)
        + Integral(in, phiHat*(x-phi)/L, quad2) ;


      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(out, phiHat*phi/L, quad2);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, phiHat, phi, vecType);

      ParameterList params;
      ParameterList linSolverParams;
      linSolverParams.set("Type", "TSF");
      linSolverParams.set("Method", "BICGSTAB");
      linSolverParams.set("Max Iterations", 1000);
      linSolverParams.set("Tolerance", 1.0e-14);
      linSolverParams.set("Precond", "ILUK");
      linSolverParams.set("Graph Fill", 1);
      linSolverParams.set("Verbosity", 4);
      
      params.set("Linear Solver", linSolverParams);
      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(params);

      Expr soln = prob.solve(linSolver);


      
      /* Project the velocity onto a discrete space so we can visualize it */
      DiscreteSpace discreteSpace(mesh, 
                                  List(new Lagrange(1), 
                                       new Lagrange(1), 
                                       new Lagrange(1)),
                                  vecType);
      L2Projector projector(discreteSpace, grad*soln);
      Expr velocity = projector.project();

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Post3d");
      w.addMesh(mesh);
      w.addField("phi", new ExprFieldWrapper(soln[0]));
      w.addField("ux", new ExprFieldWrapper(velocity[0]));
      w.addField("uy", new ExprFieldWrapper(velocity[1]));
      w.addField("uz", new ExprFieldWrapper(velocity[2]));
      w.write();


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
