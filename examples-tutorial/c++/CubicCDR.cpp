#include "Sundance.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvaluator.hpp"

using SundanceCore::List;

/** 
 * Solves a steady-state convection-diffusion equation with a cubic
 * reaction term.
 */



/* These macros define predicate classes to be used in identifying
 * boundary regions */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


static Time& sensTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("sensitivity solve"); 
  return *rtn;
}

static Time& outputTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("output"); 
  return *rtn;
}
static Time& xmlTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("xml reading"); 
  return *rtn;
}

static Time& paramsTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("xml to params conversion"); 
  return *rtn;
}

static Time& linSolveBuilderTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("building linear solver"); 
  return *rtn;
}

int main(int argc, void** argv)
{
  
  try
		{
      /* =================================================================== 
       *   Initialization  
       * =================================================================== */

      /* We will read input parameters from an XML file. The name of
       * the file is read from the command line:
       * --input=<filename>. The default is CDR.xml. */
      string filename = "CDR.xml";
      string path = "../../../examples-tutorial/c++/";
      Sundance::setOption("input", filename, "name of XML input file");
      Sundance::setOption("path", path, "path to file");
      
      /* Initialize everything, reading command-line arguments */
      Sundance::init(&argc, &argv);
      
      /* Get the processor rank. We will only write output if we 
       * are on processor 0 */
      int myRank = MPIComm::world().getRank();

      /* =================================================================== 
       *   Parameter handling 
       *
       * Here we read parameters from the XML input file, using them to build 
       * the mesh, create solver parameters, and define model parameters.
       *
       * =================================================================== */

      /* Read the parameters from the input file given on the command line */
      //      if (myRank==0) cout << "reading input file..." << endl;
      FileInputSource fis(path + filename);
      XMLObject xml;
      {
        TimeMonitor xtimer(xmlTimer());
        xml = fis.getObject();
      }

      XMLParameterListReader paramsReader;
      ParameterList params;
      {
        TimeMonitor ptimer(paramsTimer());
        params = paramsReader.toParameterList(xml);
      }

      /* Read solver parameters. The solver object will be built later, after
       * the discrete problem has been created.  */
      ParameterList solverParams = params.sublist("Nonlinear Solver");
      ParameterList linSolverParams = params.sublist("Sensitivity Solver");
      LinearSolver<double> linSolver;

      {
        TimeMonitor lstimer(linSolveBuilderTimer());
        linSolver = LinearSolverBuilder::createSolver(linSolverParams);
      }


      /* Create the mesh as specified by the Mesh parameter section */
      //      if (myRank==0) cout << "getting mesh..." << endl;

      ParameterList meshParams = params.sublist("Mesh");
      Mesh mesh = MeshBuilder::createMesh(meshParams);
      int dimension = mesh.spatialDim();

      /* Create the model parameters */
      ParameterList modelParams = params.sublist("Model");

      double alpha0 = getParameter<double>(modelParams, "Alpha");
      double beta0 = getParameter<double>(modelParams, "Beta");
      double D = getParameter<double>(modelParams, "Diffusivity");
      double L0 = getParameter<double>(modelParams, "BC Wavelength");
      double u0 = getParameter<double>(modelParams, "Peak Velocity");
      int bcQuadOrder = getParameter<int>(modelParams, "BC Quadrature Order");

      double gammaBC = 1.0;
      bool nodalBC = getParameter<bool>(modelParams, "Apply BCs at Nodes");
      bool robinBC = getParameter<bool>(modelParams, "Use Robin BCs");
      if (robinBC)
        {
          gammaBC = getParameter<double>(modelParams, "Robin BC gamma");
        }
      
      /* Read the filename for viz output */
      string outfile = getParameter<string>(params, "Output Filename");
      bool doSensitivities 
        = getParameter<bool>(params, "Compute Sensitivities");
      bool writeEachContinuationStep 
        = getParameter<bool>(params, "Write Each Continuation Step");


      /* =================================================================== 
       *   Equation definition 
       *
       * Here we write the symbolic equations and boundary conditions
       *
       * =================================================================== */

      //      if (myRank==0) cout << "setting up equation..." << endl;

      const double pi = 4.0*atan(1.0);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);


      /* Create cell filters */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter nodes = new DimensionalCellFilter(0);

      CellFilter left;
      if (nodalBC) left = nodes.subset(new LeftPointTest());
      else left = nodes.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      BasisFamily L1 = new Lagrange(1);
      Expr psi = new UnknownFunction(L1);
      Expr vPsi = new TestFunction(L1);

      /* Create gradient operator for the problem dimension */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad;
      if (dimension==1) grad = dx;
      else if (dimension==2) grad = List(dx, dy);
      else grad = List(dx, dy, dz);

      /* Set up quadrature rule. The considerations for choosing rules are:
       *
       * (*) Fourth-order quadrature is necessary to handle the cubic
       * term times the test function. We will use this on the
       * interior. While it is possible to use lower-order quadrature
       * for the other interior terms, it is most efficient to use the
       * same quadrature order for all terms; this is because using
       * the same quadrature rule lets us reuse the evaluation of the
       * discretized psi in all terms. 
       *
       * (*) Reference-integrable terms
       * are identified automatically, and the minimal-order
       * quadrature rule needed for exact integration is used in
       * computing integrals on reference elements, regardless of the
       * rule passed to the Integral() operator.  
       *
       * (*) The appropriate
       * rule for the BCs depends on the form of the velocity field
       * (for the outflow BCs) and the value of psi0 (for the
       * Dirichlet BCs).  Therefore, we read the quadrature order for
       * the BCs as an input parameter.
       */
      QuadratureFamily quad = new GaussianQuadrature(4);
      QuadratureFamily bcQuad = new GaussianQuadrature(bcQuadOrder);

      /* Define velocity field. Keep it simple for now... */
      Expr ux = u0*y*(1.0-y);
      Expr uy = 0.0;
      Expr uz = 0.0;
      Expr u;
      if (dimension==1) u = ux;
      else if (dimension==2) u = List(ux, uy);
      else u = List(ux, uy, uz);

      /* Create parameter expressions for alpha and beta. These will be
       * modified during the continuation loop, and thus must be
       * parameter exprs rather than constants. */
      Expr alpha = new Parameter(alpha0);
      Expr beta = new Parameter(beta0);

      /* Define the Ginzburg-Landau equations with an advection term added.
       * We have integrated the advection term by parts. In its present for it
       * assumes implicitly that div(u)=0, i.e., that the flow is 
       * incompressible.
       * 
       * For the outflow BC we use d_psi/dn=0, which after integration
       * by part of the advection term gives us a surface integral
       * vPsi*psi*(u_normal) on the outflow surface.
       */
      Expr eqn = Integral(interior, (grad*vPsi)*(D*(grad*psi) - u*psi), quad)
        + Integral(interior, vPsi * psi * (alpha + beta * psi*psi), quad)
        + Integral(right, vPsi*ux*psi, bcQuad);

      Expr h = new CellDiameterExpr();
      if (robinBC) eqn = eqn + Integral(left, (1.0/h)*gammaBC*vPsi*(psi - sin(2.0*pi*y/L0)), bcQuad);
        
      /* Define the Dirichlet BC */
      Expr bc;
      Expr bcScale;
      if (!robinBC)
        {
          if (nodalBC) bcScale = 1.0;
          else bcScale = 1.0/h;
          bc = EssentialBC(left, bcScale*vPsi*(psi - sin(2.0*pi*y/L0)), bcQuad);
        }



      /* ================================================================== 
       *   Discrete problem setup
       *
       * Here we create a discrete nonlinear problem
       *
       * ================================================================== */

      //      if (myRank==0) cout << "building discrete problem..." << endl;

      /* Specify that we will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create an expression for thie initial guess */
      DiscreteSpace discSpace(mesh, L1, vecType);
      Expr psi0 = new DiscreteFunction(discSpace, 0.0, "psi0");


      /* ================================================================== 
       *   Sensitivity problem setup
       *
       * Here we create problems for the sensitivities wrt alpha and beta.
       * 
       * ================================================================== */
      
      Expr w = gammaBC/h;
      Expr alphaSensEqn = Integral(interior, 
                                   (grad*vPsi)*(D*(grad*psi) - u*psi), quad)
        + Integral(interior, vPsi*(psi0 + alpha*psi 
                                   + 3.0*beta*psi0*psi0*psi), quad);
      if (robinBC) alphaSensEqn = alphaSensEqn + Integral(left, w*vPsi*psi, bcQuad);

      Expr alphaSensBC;
      if (!robinBC)
          alphaSensBC = EssentialBC(left, bcScale*vPsi*psi, bcQuad);
      
      Expr betaSensEqn = Integral(interior, 
                                  (grad*vPsi)*(D*(grad*psi) - u*psi), quad)
        + Integral(interior, vPsi*(alpha*psi0 + psi0*psi0*psi0
                                   + 3.0*beta*psi0*psi0*psi), quad);
      if (robinBC) betaSensEqn = betaSensEqn + Integral(left, w*vPsi*psi, bcQuad);
      Expr betaSensBC;
      if (!robinBC)
          betaSensBC = EssentialBC(left, bcScale*vPsi*psi, bcQuad);

      LinearProblem alphaSensProb(mesh, alphaSensEqn, alphaSensBC,
                                  vPsi, psi, vecType);

      LinearProblem betaSensProb(mesh, betaSensEqn, betaSensBC,
                                 vPsi, psi, vecType);
                   

      
      

      
      /* Create a (TSF) NonlinearOperator object representing the function F in
      * the nonlinear equation F(psi)=0 */
      NonlinearOperator<double> F 
        = new NonlinearProblem(mesh, eqn, bc, vPsi, psi, psi0, vecType);


      /* =================================================================== 
       *   Solve the problem
       *
       * We create a NOX solver using the "Nonlinear Solver" parameters read
       * in from XML. 
       *
       * The solution will overwrite the initial guess expression psi0. 
       *
       * =================================================================== */

      //      if (myRank==0) cout << "solving nonlinear system..." << endl;

      NOXSolver solver(solverParams, F);
      int numContinuationSteps 
        = getParameter<int>(solverParams, 
                            "Number of Continuation Steps");


      for (int i=1; i<=numContinuationSteps; i++)
        {
          //          if (myRank==0) cout << "continuation step " << i << endl;
          double a = i*alpha0/((double) numContinuationSteps);
          double b = i*beta0/((double) numContinuationSteps);
          alpha.setParameterValue(a);
          beta.setParameterValue(b);

          NOX::StatusTest::StatusType status = solver.solve();

          Expr sa0;
          Expr sb0;
          
          if (doSensitivities)
            {
              TimeMonitor stimer(sensTimer());
              sa0 = alphaSensProb.solve(linSolver);
              sb0 = betaSensProb.solve(linSolver);
            }
          
          if (writeEachContinuationStep)
            {
              TimeMonitor otimer(outputTimer());
              //              if (myRank==0) cout << "writing VTK output..." << endl;
              FieldWriter w = new VTKWriter(outfile + "-" 
                                            + Teuchos::toString(i));
              w.addMesh(mesh);
              w.addField("psi0", new ExprFieldWrapper(psi0[0]));
              if (doSensitivities)
                {
                  w.addField("sens_alpha", new ExprFieldWrapper(sa0[0]));
                  w.addField("sens_beta", new ExprFieldWrapper(sb0[0]));
                }
              w.write();
            }
        }




      /* =================================================================== 
       *   Viz output
       *
       * Write the solution to a VTK file. Parallelism is handled automatically. 
       *
       * =================================================================== */
      
      if (!writeEachContinuationStep)
        {
          TimeMonitor otimer(outputTimer());
          //          if (myRank==0) cout << "writing VTK output..." << endl;
          
          FieldWriter w = new VTKWriter(outfile);
          w.addMesh(mesh);
          w.addField("psi0", new ExprFieldWrapper(psi0[0]));
          w.write();
        }

      /* all done! */
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
}
