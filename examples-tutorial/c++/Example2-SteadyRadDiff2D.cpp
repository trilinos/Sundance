/* @HEADER@ */
/* @HEADER@ */

#include "Sundance.hpp"

/* ------------------------------------------------------------------------ 
 *
 * This program solves a nonlinear boundary-value problem in two spatial
 * dimensions using Galerkin finite elements and Newton's method
 *
 * ------------------------------------------------------------------------ */

const double pi = 4.0*atan(1.0);



class WestEdgeTest : public CellPredicateFunctorBase, 
                     public Playa::Handleable<CellPredicateFunctorBase>
{
public:
  WestEdgeTest() : CellPredicateFunctorBase("WestEdgeTest") {}
  virtual ~WestEdgeTest() {}
  virtual bool operator()(const Point& x) const {return fabs(x[0])<1.0e-10;}
  GET_RCP(CellPredicateFunctorBase);
};


CELL_PREDICATE(EastEdgeTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(SouthEdgeTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(NorthEdgeTest, {return fabs(x[1]-1.0) < 1.0e-10;})

bool inlineNewtonSolve(NonlinearProblem prob,
  Expr uNewt,
  int maxIters, 
  double newtTol);

bool noxNewtonSolve(const NonlinearProblem& prob, Expr uNewt);

bool playaNewtonArmijoSolve(const NonlinearProblem& prob, Expr uNewt);

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 128;
      MeshSource mesher = new PartitionedRectangleMesher(
        0.0, 1.0, nx, 1, 
        0.0, 1.0, nx, 1, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter northEdge = edges.subset(new NorthEdgeTest());
      CellFilter southEdge = edges.subset(new SouthEdgeTest());
      CellFilter eastEdge = edges.subset(new EastEdgeTest());
      CellFilter westEdge = edges.subset(new WestEdgeTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      BasisFamily bas = new Lagrange(1);
      Expr u = new UnknownFunction(bas, "u");
      Expr v = new TestFunction(bas, "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* The initial guess is u(x,y)=1 */
      DiscreteSpace discSpace(mesh, bas, vecType);
      Expr uNewt = new DiscreteFunction(discSpace, 1.0, "uNewt");

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      /* Define the weak form */
      Expr eqn = Integral(interior, u*u*u*(grad*v)*(grad*u), quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(northEdge, v*(u - pow(1.0 + sin(pi*x),0.25)),quad)
        + EssentialBC(southEdge + eastEdge + westEdge, 
          v*(u - 1.0),quad);

      /* We can now set up the nonlinear problem! */
      NonlinearProblem prob(mesh, eqn, bc, v, u, uNewt, vecType); 

      bool nonlinSolveOK = false;
      bool useNox = false;
      if (useNox)
      {
        nonlinSolveOK = noxNewtonSolve(prob, uNewt);
      }
      else
      {
        nonlinSolveOK = playaNewtonArmijoSolve(prob, uNewt);
      }
      
      TEUCHOS_TEST_FOR_EXCEPT(!nonlinSolveOK);

      FieldWriter writer = new VTKWriter("SteadyRadDiff2D"); 
      writer.addMesh(mesh);
      writer.addField("u", new ExprFieldWrapper(uNewt[0]));
      writer.write();

      Expr uExact = pow(1.0 + sin(pi*x)*sinh(pi*y)/sinh(pi), 0.25);
      double err = L2Norm(mesh, interior, uNewt-uExact, quad);
      Out::root() << "error = " << setw(16) << err << endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(err, tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}


bool inlineNewtonSolve(NonlinearProblem prob,
  Expr uNewt,
  int maxIters, 
  double newtTol)
{
  bool newtonConverged = false;
  LinearSolver<double> linSolver 
    = LinearSolverBuilder::createSolver("amesos.xml");

  /* Allocate objects for the Jacobian, residual, and Newton step */
  LinearOperator<double> J = prob.allocateJacobian();
  Vector<double> resid = J.range().createMember();
  Vector<double> newtonStep = J.domain().createMember();

  for (int i=0; i<maxIters; i++)
  {
    prob.setEvalPoint(uNewt);
    prob.computeJacobianAndFunction(J, resid);
    SolverState<double> solveState 
      = linSolver.solve(J, -1.0*resid, newtonStep);
    
    TEUCHOS_TEST_FOR_EXCEPTION(solveState.finalState() != SolveConverged,
      std::runtime_error,
      "linear solve failed!");
    
    addVecToDiscreteFunction(uNewt, newtonStep);
    double newtStepNorm = newtonStep.norm2();
    Out::root() << "|newt step| = " << newtStepNorm << endl;
    if (newtStepNorm < newtTol) 
    {
      newtonConverged = true;
      break;
    }
  }

  return newtonConverged;
}

bool noxNewtonSolve(const NonlinearProblem& prob, Expr uNewt)
{
  /* Use a NOX nonlinear solver */
  ParameterXMLFileReader reader("nox-amesos.xml");
  ParameterList solverParams = reader.getParameters();
  NOXSolver solver(solverParams);
  
  SolverState<double> stat = prob.solve(solver);

  return stat.finalState()==SolveConverged;
}

bool playaNewtonArmijoSolve(const NonlinearProblem& prob, Expr uNewt)
{
  /* Use the Playa Newton-Armijo nonlinear solver */
  LinearSolver<double> linSolver = LinearSolverBuilder::createSolver("amesos.xml");
  ParameterList solverParams("NewtonArmijoSolver");
  solverParams.set("Tau Relative", 1.0e-12);
  solverParams.set("Tau Absolute", 1.0e-12);
  solverParams.set("Alpha", 1.0e-4);
  solverParams.set("Verbosity", 3);

  NonlinearSolver<double> solver = new NewtonArmijoSolver<double>(solverParams, linSolver);
  
  SolverState<double> stat = prob.solve(solver);
  if (stat.finalState() != SolveConverged)
  {
    Out::os() << stat.finalMsg() << endl;
  }

  return stat.finalState()==SolveConverged;
}
