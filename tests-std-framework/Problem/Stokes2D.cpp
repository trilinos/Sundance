#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

/** 
 * Solves the Stokes equation in 2D
 */

bool leftPointTest(const Point& x) {return fabs(x[0]+1.0) < 1.0e-10;}
bool bottomPointTest(const Point& x) {return fabs(x[1]+1.0) < 1.0e-10;}
bool rightPointTest(const Point& x) {return fabs(x[0]-1.0) < 1.0e-10;}
bool topPointTest(const Point& x) {return fabs(x[1]-1.0) < 1.0e-10;}

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
      int nx = 40;
      int ny = 40;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(-1.0, 1.0, nx*np, np,
                                                         -1.0, 1.0, ny, 1,
                                                         meshType);




      Mesh mesh = mesher.getMesh();
      double h = 1.0/((double) ny);

//       FieldWriter wMesh = new VerboseFieldWriter();
//       wMesh.addMesh(mesh);
//       wMesh.write();

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
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(1), "u_y");
      Expr vy = new TestFunction(new Lagrange(1), "v_y");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

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
      double beta = 0.2;
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy) - p*(dx*vx+dy*vy)
                          + h*h*beta*(grad*q)*(grad*p) - q*(dx*ux+dy*uy),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = 0.5*(1.0-y*y);
      Expr bc = EssentialBC(left, vx*ux + vy*uy , quad4)
        + EssentialBC(top, vx*(ux-y) + vy*uy, quad2)
        + EssentialBC(right, vx*ux + vy*uy, quad2)
        + EssentialBC(bottom, vx*ux + vy*uy, quad2);


    //   Expr poissonEqn = Integral(interior, (grad*vx)*(grad*ux), quad2);
//       Expr poissonBC = EssentialBC(left, vx*(ux-uInflow), quad2)
//         + EssentialBC(top, vx*ux, quad2)
//         + EssentialBC(bottom, vx*ux, quad2);

      
//       Expr ppEqn = Integral(interior, h*h*beta*(grad*vx)*(grad*ux), quad2);
//       Expr ppBC;

//       Expr xConstraint = Integral(interior, vx*dx*ux, quad2);
//       Expr xConstraintBC;
//       Expr yConstraint = Integral(interior, vx*dy*ux, quad2);
//       Expr yConstraintBC;

//       Expr xP = Integral(interior, ux*dx*vx, quad2);
//       Expr xPBC;
//       Expr yP = Integral(interior, ux*dy*vx, quad2);
//       Expr yPBC;

//       LinearProblem poissonProb(mesh, poissonEqn, poissonBC, vx, ux, vecType);
//       LinearProblem ppProb(mesh, ppEqn, ppBC, vx, ux, vecType);
//       LinearProblem xProb(mesh, xConstraint, xConstraintBC, vx, ux, vecType);
//       LinearProblem yProb(mesh, yConstraint, yConstraintBC, vx, ux, vecType);
//       LinearProblem xpProb(mesh, xP, xPBC, vx, ux, vecType);
//       LinearProblem ypProb(mesh, yP, yPBC, vx, ux, vecType);

      Assembler::workSetSize() = 100;
      FunctionalEvaluator::workSetSize() = 100;
      //      Assembler::classVerbosity() = VerbExtreme;

 //      cerr << "--------------- Poisson operator " << endl;
//       LinearOperator<double> A_poisson = poissonProb.getOperator();
//       A_poisson.print(cerr);
//       cerr << "--------------- Poisson vector " << endl;
//       Vector<double> b_poisson = poissonProb.getRHS();
//       cerr << b_poisson << endl;

      

//       cerr << "---------------- pressure Poisson operator " << endl;
//       LinearOperator<double> A_pp = ppProb.getOperator();
//       A_pp.print(cerr);

//       cerr << "--------------- x constraint operator " << endl;
//       LinearOperator<double> A_x = xProb.getOperator();
//       A_x.print(cerr);


//       cerr << "--------------- y constraint operator " << endl;
//       LinearOperator<double> A_y = yProb.getOperator();
//       A_y.print(cerr);

//       cerr << "--------------- x pressure force operator " << endl;
//       LinearOperator<double> A_px = xpProb.getOperator();
//       A_px.print(cerr);


//       cerr << "--------------- y pressure force operator " << endl;
//       LinearOperator<double> A_py = ypProb.getOperator();
//       A_py.print(cerr);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q), 
                         List(ux, uy, p), vecType);

  //     cerr << "------------------ stokes operator " << endl;
//       LinearOperator<double> A_stokes = prob.getOperator();
//       A_stokes.print(cerr);
      
//       cerr << "--------------- Stokes vector " << endl;
//       Vector<double> b_stokes = prob.getRHS();
//       cerr << b_stokes << endl;


      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 2;
      azParams[AZ_max_iter] = 1000;
      azParams[AZ_tol] = 1.0e-10;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Stokes2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
