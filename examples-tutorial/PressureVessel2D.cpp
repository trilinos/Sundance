/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 2D
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
        = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/vessel2D.ncdf", meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter inside = edges.labeledSubset(2);
      CellFilter outside = edges.labeledSubset(1);
      CellFilter xNormalFace = edges.labeledSubset(3);
      CellFilter yNormalFace = edges.labeledSubset(4);

      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
			Expr vx = new TestFunction(new Lagrange(2), "vx");
			Expr vy = new TestFunction(new Lagrange(2), "vy");
			Expr ux = new UnknownFunction(new Lagrange(2), "ux");
			Expr uy = new UnknownFunction(new Lagrange(2), "uy");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      
      /* Young's modulus */
      double E = 1.0;
      /* Poisson's ratio */
      double nu = 0.25;

      double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
      double mu = E/2.0/(1.0 + nu);

      /* form strain tensors using voigt notation */
			Expr strain = List(dx*ux, dy*uy, dx*uy + dy*ux);
			Expr varStrain = List(dx*vx, dy*vy, dx*vy + dy*vx);

			/* */
			Expr D = List(List(lambda+2.0*mu,    lambda,             0.0), 
                    List(lambda,           lambda+2.0*mu,      0.0),
                    List(0.0,              0.0,                mu));

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      /* Define the weak form */
      Expr eqn = Integral(interior, varStrain*(D*strain), quad2)
        + Integral(inside, -x*vx, quad2);

      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(xNormalFace, vx*ux, quad2)
        + EssentialBC(yNormalFace, vy*uy, quad2);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy), 
                         List(ux, uy), vecType);

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azOptions[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-10;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Elasticity2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.write();


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
