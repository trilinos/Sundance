#include "Sundance.hpp"

/* We solve in 2D
   - Laplace u = 0
   stationary
   With Dirichlet boundary condition:
             u = const , at the boundary
 */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-4;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-4;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-4;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-4;})

int main(int argc, char** argv)
{
  try
  {
      Sundance::init(&argc, &argv);
      int myrank = MPIComm::world().getRank();

      // We will do our linear algebra using Epetra 
      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new HNodeMeshType2D();
      MeshSource mesher = new HNodeMesher2D(0.0, 0.0, 1.0 , 1.0 , 0.1 , 0.1, meshType);
      Mesh mesh = mesher.getMesh();
   
      // set the Verbosity high for the grid
      //Mesh::classVerbosity() = VerbExtreme; 
    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 0);
    watchMe.setParam("integration setup", 1);
    watchMe.setParam("integral transformation", 1);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 0);
    watchMe.setParam("evaluation", 0);


      cout << "Nr Points  "<<mesh.numCells(0) << endl;
      cout << "My Rank is :" << myrank << endl;
      
      // Define the domains 
      CellFilter interior = new MaximalCellFilter();
      CellFilter boundary = new BoundaryCellFilter();

      // Find left, right, top, bottom boundary edges
      CellFilter left = boundary.subset(new LeftPointTest());
      CellFilter right = boundary.subset(new RightPointTest());
      CellFilter top = boundary.subset(new TopPointTest());
      CellFilter bottom = boundary.subset(new BottomPointTest());

      BasisFamily La1 = new Lagrange(1);

      // Define the unknown function and the testfunction 
      Expr vs = new TestFunction( La1 , "v" );
      Expr us = new UnknownFunction( La1 , "u" );

      // Define the derivatives
      Expr dx = new Derivative(0); 
      Expr dy = new Derivative(1);
      Expr grad = List(dx,dy);

      // Define the quadrature
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      // Define the integral 
      Expr eqn = Integral(interior, (grad*us)*(grad*vs) - vs*3.0, quad2 );

      // Define the boundary conditions
      Expr bc =  EssentialBC( right , vs*(us-1.0) , quad2 );

      // We can now set up the linear problem!
      LinearProblem prob(mesh, eqn, bc, vs, us, vecType);

      // Read the parameters for the linear solver from an XML file

      ParameterXMLFileReader reader("bicgstab.xml");
      //ParameterXMLFileReader reader("superlu_dist.xml");

      ParameterList solverParams = reader.getParameters();
      // Now we define the problem
      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(solverParams);
 
      // solve the problem
      Expr soln = prob.solve(linSolver);
      // Write the field in VTK format 
      FieldWriter w = new VTKWriter("Poisson_HN");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(soln[0]));
      w.write(); 
      
      /* ---------- TESTING SOLUTION -------- */

      Expr checkLeftBoundaryErr = Integral( left, (soln[0]-2.5)*(soln[0]-2.5),
       		new GaussianQuadrature(2) );

      FunctionalEvaluator errDif(mesh, checkLeftBoundaryErr);

      double leftBoundaryerror = errDif.evaluate();

      std::cout << "leftBoundaryerror = " << leftBoundaryerror << std::endl;

      double tol = 5.0e-2;
      Sundance::passFailTest(sqrt(leftBoundaryerror), tol);

  }
  catch(exception& e)
  {
      Sundance::handleException(e);
  }
  Sundance::finalize();
}
