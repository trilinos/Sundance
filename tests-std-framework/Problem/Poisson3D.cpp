#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"

/** 
 * Solves the Poisson equation in 2D
 */


int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      VerbositySetting verb = VerbSilent;
      Assembler::classVerbosity() = verb;

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/cube.ncdf", meshType);
      Mesh mesh = mesher.getMesh();

      if (verb > VerbHigh)
        {
          FieldWriter wMesh = new VerboseFieldWriter();
          wMesh.addMesh(mesh);
          wMesh.write();
        }

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);
      CellPredicate topFunc = new LabelCellPredicate(1);
      CellPredicate bottomFunc = new LabelCellPredicate(2);
      CellFilter top = faces.subset(topFunc);
      CellFilter bottom = faces.subset(bottomFunc);

      cerr << "top cells = " << top.getCells(mesh) << endl;
      cerr << "bottom cells = " << bottom.getCells(mesh) << endl;

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
      Expr eqn = Integral(interior, (grad*v)*(grad*u)  + v, quad2);

      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, v*u, quad2)
        + EssentialBC(top, v*(u-z), quad2);

      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;

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

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Poisson3d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.write();


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
