#include "Sundance.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;
/** 
 * \example PressureStabilizedNS2D.cpp
 * 
 * Solves the Navier-Stokes equation in 2D using pressure stabilization.
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
      int nx = 64;
      int ny = 64;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx*np, np,
                                                         0.0, 1.0, ny, 1,
                                                         meshType);




      Mesh mesh = mesher.getMesh();
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      //    Expr h = 2.0/((double) nx); 
      //Expr h = new CellDiameterExpr();
      double h0 = 2.0/((double) nx); 
      // Expr h = h0*pow(1.0+x, 0.0);
      Expr h = h0;

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter points = new DimensionalCellFilter(0);
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
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(1), "u_y");
      Expr vy = new TestFunction(new Lagrange(1), "v_y");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");
      Expr u = List(ux, uy);

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad1 = new GaussianQuadrature(1);
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Stabilization parameters */
      double beta = 4.0 * 0.02;

      /* A parameter expression for the Reynolds number */
      Expr reynolds = new Parameter(20.0);

      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + beta*h*h*(grad*q)*(grad*p) + q*(dx*ux+dy*uy),
                          quad1)
        + Integral(interior, reynolds*(vx*(u*grad)*ux)
                   + reynolds*(vy*(u*grad)*uy), quad2);
        
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(left, vx*ux + vy*uy, quad2)
        + EssentialBC(right, vx*ux + vy*uy, quad2)
        + EssentialBC(top, vx*(ux-1.0) + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, SundanceStdFwk::List(L1, L1, L1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");
      
      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, List(vx, vy, q),
                               List(ux, uy, p), u0, vecType);

      

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
          FieldWriter w = new VTKWriter("ns-r" + Teuchos::toString(Re));
          w.addMesh(mesh);
          w.addField("u_x", new ExprFieldWrapper(u0[0]));
          w.addField("u_y", new ExprFieldWrapper(u0[1]));
          w.addField("p", new ExprFieldWrapper(u0[2]));
          w.write();
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
