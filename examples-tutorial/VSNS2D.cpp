#include "Sundance.hpp"

/** 
 * Solves the Navier-Stokes equations on the lid-driven cavity
 */

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Read the mesh */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../examples-tutorial/square128.ncdf", meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter bottom = edges.labeledSubset(1);
      CellFilter right = edges.labeledSubset(2);
      CellFilter top = edges.labeledSubset(3);
      CellFilter left = edges.labeledSubset(4);

      
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

      /* A parameter expression for the Reynolds number */
      Expr reynolds = new Parameter(1.0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad1 = new GaussianQuadrature(1);
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      /* Define the weak form */
      Expr psiEqn = Integral(interior, (grad*vPsi)*(grad*psi) + vPsi*omega, 
                             quad2)
        + Integral(top, -1.0*vPsi, quad1);

      Expr omegaEqn = Integral(interior, (grad*vOmega)*(grad*omega)
                               + reynolds*vOmega*((dy*psi)*dx*omega 
                                                  - (dx*psi)*dy*omega),
                               quad1);
      
      Expr eqn = omegaEqn + psiEqn;

      /* Define the Dirichlet BC */
      CellFilter walls = left + bottom + right;
      Expr bc = EssentialBC(walls, vOmega*psi, quad2) 
        + EssentialBC(top, vOmega*psi, quad2) ;

      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, List(L1, L1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, List(vPsi, vOmega),
                               List(psi, omega), u0, vecType);

      ParameterXMLFileReader reader("../../examples-tutorial/nox.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams, F);

      int numReynolds = 10;
      double finalReynolds = 500.0;
      for (int r=1; r<=numReynolds; r++)
        {
          double Re = r*finalReynolds/((double) numReynolds);
          reynolds.setParameterValue(Re);
          cerr << "--------------------------------------------------------- " << endl;
          cerr << " solving for Reynolds Number = " << Re << endl;
          cerr << " reynolds = " << reynolds << endl;
          cerr << "--------------------------------------------------------- " << endl;
          // Solve the nonlinear system
          NOX::StatusTest::StatusType status = solver.solve();

          /* Write the field in VTK format */
          FieldWriter w = new VTKWriter("vns-r" + Teuchos::toString(Re));
          w.addMesh(mesh);
          w.addField("streamfunction", new ExprFieldWrapper(u0[0]));
          w.addField("vorticity", new ExprFieldWrapper(u0[1]));
          w.write();
        }
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
