#include "Sundance.hpp"
#include "PDEOptNonlinearPDEConstrainedObj.hpp"

#include "PlayaBasicLMBFGS.hpp"
#include "PlayaOptBuilder.hpp"


int main(int argc, char** argv)
{
  try
		{
			Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      
      int nx = 2;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();
      CellFilter interior = new MaximalCellFilter();
      
      /* Create a vector space factory, used to 
       * specify the low-level linear algebra representation */
      VectorType<double> vecType = new EpetraVectorType();
  
      /* create a discrete space on the mesh */
      DiscreteSpace discreteSpace(mesh, new Lagrange(1), vecType);

      /* initialize the design, state, and multiplier vectors */
      Expr alpha0 = new DiscreteFunction(discreteSpace, 0.5, "alpha0");
      Expr u0 = new DiscreteFunction(discreteSpace, 0.5, "u0");
      Expr lambda0 = new DiscreteFunction(discreteSpace, 0.5, "lambda0");

      /* create symbolic objects for test and unknown functions */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr lambda = new UnknownFunction(new Lagrange(1), "lambda");
      Expr alpha = new UnknownFunction(new Lagrange(1), "alpha");

      /* create quadrature rules of different orders */
      QuadratureFamily q1 = new GaussianQuadrature(1);
      QuadratureFamily q2 = new GaussianQuadrature(2);
      QuadratureFamily q4 = new GaussianQuadrature(4);

      /* Regularization weight */
      double R = 1.0;

      /* Exact solution (computed by Mathematica to 6 digits) */
      double U0 = 1.20404;
      double A0 = -0.20863;
      Expr aExact = new DiscreteFunction(discreteSpace, A0, "alpha0");
      Expr uExact = new DiscreteFunction(discreteSpace, U0, "u0");

      /* Form objective function */
      Expr reg = Integral(interior, 0.5 * R * alpha*alpha, q2);
      Expr fit = Integral(interior, 0.5 * pow(u-1.0, 2.0), q2);

      Expr constraintEqn = Integral(interior, 
        lambda*(alpha - sin(u-sqrt(2.0))), q4);
      Expr L = reg + fit + constraintEqn;

      Expr constraintBC;

      Functional Lagrangian(mesh, L, constraintBC, vecType);
      
      LinearSolver<double> adjSolver 
        = LinearSolverBuilder::createSolver("amesos.xml");
      ParameterXMLFileReader reader("nox.xml");
      ParameterList noxParams = reader.getParameters();
      NOXSolver nonlinSolver(noxParams);

      RCP<ObjectiveBase> obj = rcp(new NonlinearPDEConstrainedObj(
        Lagrangian, u, u0, lambda, lambda0, alpha, alpha0,
        nonlinSolver, adjSolver));

      Vector<double> xInit = obj->getInit();

      bool doFDCheck = true;
      if (doFDCheck)
      {
        Out::root() << "Doing FD check of gradient..." << endl;
        bool fdOK = obj->fdCheck(xInit, 1.0e-6, 2);
        if (fdOK) 
        {
          Out::root() << "FD check OK" << endl;
        }
        else
        {
          Out::root() << "FD check FAILED" << endl;
        }
      }

      RCP<UnconstrainedOptimizerBase> opt 
          = OptBuilder::createOptimizer("basicLMBFGS.xml");
      opt->setVerb(3);

      OptState state = opt->run(obj, xInit);

      bool ok = true;
      if (state.status() != Opt_Converged)
      {
        Out::root()<< "optimization failed: " << state.status() << endl;
        ok = false;
      }
      else
      {
        Out::root() << "opt converged: " << state.iter() << " iterations"
                    << endl;
        Out::root() << "exact solution: U0=" << U0 << endl;
      }
      FieldWriter w = new MatlabWriter("AlgebraicConstrained1D.dat");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0));
      w.addField("alpha", new ExprFieldWrapper(alpha0));
      w.addField("lambda", new ExprFieldWrapper(lambda0));
      w.write();

          
      double uErr = L2Norm(mesh, interior, u0-uExact, q4);
      double aErr = L2Norm(mesh, interior, alpha0-aExact, q4);
      Out::root() << "error in u = " << uErr << endl;
      Out::root() << "error in alpha = " << aErr << endl;

      double tol = 0.001;
      Sundance::passFailTest(uErr + aErr, tol);
      
      
    }
	catch(exception& e)
		{
      cerr << "main() caught exception: " << e.what() << endl;
		}
	Sundance::finalize();
  return Sundance::testStatus(); 
}
