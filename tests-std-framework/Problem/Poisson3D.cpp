#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"

using SundanceCore::List;
/** 
 * Solves the Poisson equation in 3D
 */


int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Read the mesh */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../../tests-std-framework/Problem/cube-coarse.ncdf", meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);
      CellFilter side1 = faces.labeledSubset(1);
      CellFilter side2 = faces.labeledSubset(2);
      CellFilter side3 = faces.labeledSubset(3);
      CellFilter side4 = faces.labeledSubset(4);
      CellFilter side5 = faces.labeledSubset(5);
      CellFilter side6 = faces.labeledSubset(6);

      
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
      Expr eqn = Integral(interior, (grad*v)*(grad*u) +2.0*v, quad2);

      /* Define the Dirichlet BC */
      Expr exactSoln = (x + 1.0)*x - 1.0/4.0;
      Expr bc = EssentialBC(side4, v*(u-exactSoln), quad4)
        + EssentialBC(side6, v*(u-exactSoln), quad4);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);


      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      L2Projector proj1(discSpace, exactSoln);
      L2Projector proj2(discSpace, soln-exactSoln);
      L2Projector proj3(discSpace, pow(soln-exactSoln, 2.0));
      Expr exactDisc = proj1.project();
      Expr errorDisc = proj2.project();
      Expr errorSqDisc = proj3.project();

      cerr << "writing fields" << endl;
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Poisson3d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("exact soln", new ExprFieldWrapper(exactDisc));
      w.addField("error", new ExprFieldWrapper(errorDisc));
      w.addField("errorSq", new ExprFieldWrapper(errorSqDisc));
      w.write();

      cerr << "computing error" << endl;

      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2.0),
                              new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cerr << "error norm = " << sqrt(errorSq) << endl << endl;

      double tol = 1.0e-12;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
