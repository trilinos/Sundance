#include "Sundance.hpp"
#include "SundanceEvaluator.hpp"

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"


#include "TSFNOXSolver.H"

/** 
 * Solves Kepler's equation x = u(x) + e*sin(u(x))
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
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      const double pi = 4.0*atan(1.0);
      MeshSource mesher = new PartitionedLineMesher(0.0, pi, 10*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

      /* Create coordinate function */
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

     
      /* Define the weak form */
      Expr ecc = new Parameter(0.5, "e");
      Expr eqn = Integral(interior, v*(u + ecc*sin(u) - x), quad);
      Expr bc;

      /* Create a TSF NonlinearOperator object */
      NonlinearOperator<double> F = new NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
      F.verbosity() = VerbExtreme;

      ParameterXMLFileReader reader("../../../tests-std-framework/Problem/nox.xml");
      ParameterList noxParams = reader.getParameters();

      cerr << "solver params = " << noxParams << endl;

      NOXSolver solver(noxParams, F);

      int numEcc = 10;
      double finalEcc = 0.95;
      for (int r=1; r<=numEcc; r++)
        {
          double e = r*finalEcc/((double) numEcc);
          ecc.setParameterValue(e);
          cerr << "--------------------------------------------------------- " << endl;
          cerr << " solving for eccentricity = " << ecc << endl;
          cerr << "--------------------------------------------------------- " << endl;
          // Solve the nonlinear system
          NOX::StatusTest::StatusType status = solver.solve();

          /* Write the field in matlab format */
          FieldWriter w = new MatlabWriter("kepler-e" + Teuchos::toString(e) + ".dat");
          w.addMesh(mesh);
          w.addField("true anomaly", new ExprFieldWrapper(u0[0]));
          w.write();
        }

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  MPISession::finalize();
}
