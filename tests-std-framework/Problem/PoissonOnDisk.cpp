#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D on the unit disk
 */

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();



      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Get a mesh */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource meshReader = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/disk.ncdf", meshType);

      Mesh mesh = meshReader.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

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
      Expr eqn = Integral(interior, (grad*v)*(grad*u)  + v, quad2);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bdry, v*u, quad4);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-10;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      double R = 1.0;
      Expr exactSoln = 0.25*(x*x + y*y - R*R);

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr du = L2Projector(discSpace, exactSoln-soln).project();

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("PoissonOnDisk");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("error", new ExprFieldWrapper(du));
      w.write();

      
      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2),
                              new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(sqrt(errorSq), tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
