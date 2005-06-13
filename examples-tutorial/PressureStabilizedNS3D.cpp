#include "Sundance.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
/** 
 * \example PressureStabilizedNS2D.cpp
 * 
 * Solves the Navier-Stokes equation in 2D using pressure stabilization.
 */


int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../examples-tutorial/box-0.05.ncdf", meshType);
      Mesh mesh = mesher.getMesh();



      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);
      //    Expr h = 2.0/((double) nx); 
      //Expr h = new CellDiameterExpr();
      double h0 = 2.0/((double) 20); 
      // Expr h = h0*pow(1.0+x, 0.0);
      Expr h = h0;

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(2);
      CellFilter points = new DimensionalCellFilter(0);

      CellFilter left = edges.labeledSubset(3);
      CellFilter right = edges.labeledSubset(4);
      CellFilter top = edges.labeledSubset(1);
      CellFilter bottom = edges.labeledSubset(2);
      CellFilter front = edges.labeledSubset(5);
      CellFilter back = edges.labeledSubset(6);
      CellFilter walls = left + right + front + back + bottom;
      CellFilter lid = top;
      CellFilter pinnedNode = points.labeledSubset(1);


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(1), "u_y");
      Expr vy = new TestFunction(new Lagrange(1), "v_y");
      Expr uz = new UnknownFunction(new Lagrange(1), "u_z");
      Expr vz = new TestFunction(new Lagrange(1), "v_z");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");
      Expr u = List(ux, uy, uz);
      Expr v = List(vx, vy, vz);

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad1 = new GaussianQuadrature(1);
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Stabilization parameters */
      double beta = 4.0 * 0.02;

      /* A parameter expression for the Reynolds number */
      Expr reynolds = new Parameter(150.0);

      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy) + (grad*vz)*(grad*uz) 
                          - p*(grad*v)
                          + beta*h*h*(grad*q)*(grad*p) + q*(grad*u), quad2)
        + Integral(interior, reynolds*(vx*(u*grad)*ux)
                   + reynolds*(vy*(u*grad)*uy) 
                   + reynolds*(vz*(u*grad)*uz), quad2);
        
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(walls, v*u, quad2)
        + EssentialBC(lid, v*u - vx, quad2)
        + EssentialBC(pinnedNode, q*p, quad1);


      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, SundanceStdFwk::List(L1, L1, L1, L1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");
      
      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, List(vx, vy, vz, q),
                               List(ux, uy, uz, p), u0, vecType);

      

      ParameterXMLFileReader reader("../../examples-tutorial/nox.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams, F);



      int numReynolds = 1;
      double finalReynolds = 10.0;
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
          FieldWriter w = new VTKWriter("ns-r" + Teuchos::toString(Re));
          w.addMesh(mesh);
          w.addField("u_x", new ExprFieldWrapper(u0[0]));
          w.addField("u_y", new ExprFieldWrapper(u0[1]));
          w.addField("u_z", new ExprFieldWrapper(u0[2]));
          w.addField("p", new ExprFieldWrapper(u0[3]));
          w.write();
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
