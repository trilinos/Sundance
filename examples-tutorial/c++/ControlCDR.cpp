/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceDefs.hpp"
#include "TSFConfigDefs.hpp"

#ifndef HAVE_ENABLED_MOOCHO

#include <iostream>

int main()
{
  std::cout << "moocho not present: test INACTIVE" << endl;
}

#else

#include "MoochoPack_MoochoSolver.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_LOWSFactoryBuilder.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"

#include "SundanceNLPModelEvaluator.hpp"
#include "Sundance.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

static Time& moochoTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("moocho solve"); 
  return *rtn;
}



/** 
 * Optimal target matching for steady-state convection-diffusion equation with a cubic
 * reaction term.
 */

/* These macros define predicate classes to be used in identifying
 * boundary regions */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


namespace Thyra
{
  class CDRControlModel : public SundanceNLPModelEvaluator
  {
  public:
    /** */
    CDRControlModel(const ParameterList& params,
                    const VectorType<double>& vecType);

    /** */
    Vector<double> getInitialState() const 
    {
      Vector<double> rtn = stateSpace().createMember();
      rtn.setToConstant(0.0);
      return rtn;
    }

    

    /** */
    const string& outfile() const {return outfile_;}

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const Expr& target() const {return target_;}

    /** */
    void setTargetVector(const Expr& expr);

  private:
    string outfile_;
    Mesh mesh_;
    Expr target_;
  };
}


int main(int argc, void** argv)
{
  using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;

  try
		{
      /* We will read input parameters from an XML file. The name of
       * the file is read from the command line:
       * --input=<filename>. The default is CDR.xml. */
      string path = "../../../examples-tutorial/c++/";
      string filename = "CDR_control.xml";
      Sundance::setOption("input", filename, "name of XML input file");
      Sundance::setOption("path", path, "path to XML input file");

      /* Initialize everything, reading command-line arguments */
      Sundance::init(&argc, &argv);

      /* Get the processor rank. We will only write output if we 
       * are on processor 0 */
      int myRank = MPIComm::world().getRank();

      /* Read the parameters from the input file given on the command line */
      if (myRank==0) cout << "reading input file..." << endl;
      FileInputSource fis(path + filename);
      XMLObject xml = fis.getObject();
      XMLParameterListReader paramsReader;
      ParameterList params = paramsReader.toParameterList(xml);

      if (myRank==0) cout << " CDR input parameters = " << endl << params << endl;

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* create an object for the physics model */
      RefCountPtr<Thyra::CDRControlModel> cdrModel 
        = rcp(new Thyra::CDRControlModel(params, vecType));


      RefCountPtr<Thyra::ModelEvaluator<double> > model = cdrModel;

      
      /* specify the linear solver to be used in Moocho */

      RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory
        = LOWSFactoryBuilder::createLOWSFactory(params.sublist("Linear Solver"));

      RefCountPtr<Thyra::ModelEvaluator<double> > modelWithLOWS
        = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model,
                                                                       lowsFactory));

      

      /* set up the solution of the forward problem to obtain a target */
      TEST_FOR_EXCEPT(!params.isSublist("Forward Problem"));
      ParameterList fwdParams = params.sublist("Forward Problem");
      Expr fwdSoln = cdrModel->solveForward(fwdParams);
      
      string outputFilename = "fwdSoln";
      if (fwdParams.isParameter("Output Filename"))
        {
          outputFilename = getParameter<string>(fwdParams, "Output Filename");
        }
      
      cdrModel->setTargetVector(fwdSoln);

      /* perturb the initial guess from the values used to make the target */
      Array<double> initialGuess;

      if (myRank==0)
        {
          initialGuess = cdrModel->paramArray(fwdParams, "Design Parameters");
          double perturbation = getParameter<double>(params, "Perturbation Magnitude");
          
          for (unsigned int i=0; i<initialGuess.size(); i++)
            {
              double a = 1.0-perturbation;
              double b = 2.0*perturbation;
              initialGuess[i] = initialGuess[i] * (a + b*rand()/((double)RAND_MAX));
            }
        }
      MPIContainerComm<double>::bcast(initialGuess, 0, cdrModel->mesh().comm());
      cdrModel->setInitialParameters(initialGuess);



      /* set up a Moocho NLP object */
      NLPFirstOrderThyraModelEvaluator nlp(modelWithLOWS, 0, 0);
      
      /* Create the optimization solver object */
      MoochoSolver  solver;
      solver.set_nlp( Teuchos::rcp(&nlp,false) );
      

      /* do it! */
      {
        TimeMonitor timer(moochoTimer());
        if (myRank==0) cout << "starting solve" << endl;
        const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();
      }

      /* write viz output */
      Expr uFinal = cdrModel->stateVariable();
      FieldWriter w2 = new VTKWriter(cdrModel->outfile());
      w2.addMesh(cdrModel->mesh());
      w2.addField("Target", new ExprFieldWrapper(fwdSoln[0]));
      w2.addField("Optimal Psi", new ExprFieldWrapper(uFinal[0]));
      w2.write();
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}


namespace Thyra
{
  CDRControlModel::CDRControlModel(const ParameterList& paramList,
                                   const VectorType<double>& vecType)
    : SundanceNLPModelEvaluator(vecType), outfile_(), mesh_(), target_()
  {
    int myRank = MPIComm::world().getRank();
    /* Create the mesh as specified by the Mesh parameter section */
    if (myRank==0) cout << "getting mesh..." << endl;
    
    ParameterList meshParams = paramList.sublist("Mesh");
    Mesh mesh = MeshBuilder::createMesh(meshParams);
    mesh_ = mesh;
    int dimension = mesh.spatialDim();
    
    /* Create the model parameters */
    ParameterList modelParams = paramList.sublist("Model");
    
    double alpha0 = getParameter<double>(modelParams, "Alpha");
    double beta0 = getParameter<double>(modelParams, "Beta");
    double D = getParameter<double>(modelParams, "Diffusivity");
    double u0 = getParameter<double>(modelParams, "Peak Velocity");
    int bcQuadOrder = getParameter<int>(modelParams, "BC Quadrature Order");
    int numControls = getParameter<int>(modelParams, "Num Controls");

    double gammaBC = 1.0;
    bool nodalBC = getParameter<bool>(modelParams, "Apply BCs at Nodes");
    bool robinBC = getParameter<bool>(modelParams, "Use Robin BCs");
    if (robinBC)
      {
        gammaBC = getParameter<double>(modelParams, "Robin BC gamma");
      }

    /* register alpha and beta as continuation parameters */
    Expr ab = List(new Parameter(alpha0), new Parameter(beta0));
    Expr abFinal = List(new Parameter(alpha0), new Parameter(beta0));
    setContinuationParameters(ab);
    setFinalContinuationValues(abFinal);
    Expr alpha = ab[0];
    Expr beta = ab[1];
    
    /* Read the filename for viz output */
    outfile_ = getParameter<string>(paramList, "Output Filename");


    /* =================================================================== 
     *   Equation definition 
     *
     * Here we write the symbolic equations and boundary conditions
     *
     * =================================================================== */
    
    if (myRank==0) cout << "setting up equation..." << endl;
    
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
    Expr dx = new SundanceCore::Derivative(0);
    Expr dy = new SundanceCore::Derivative(1);
    Expr dz = new SundanceCore::Derivative(2);
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
    
    
    /* Define the boundary field in terms of the design parameters */
    Expr control = 0.0;
    Array<Expr> p(numControls);
    for (unsigned int i=1; i<=p.size(); i++) 
      {
        p[i-1] = new Parameter(1.0);
        control = control + p[i-1]*sin(i*pi*y);
      }
    Expr param = new ListExpr(p);

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
        
    /* Define the left BC using the control field */
    Expr bc;
    Expr bcScale;
    Expr h = new CellDiameterExpr();
    if (!robinBC)
      {
        if (nodalBC) bcScale = 1.0;
        else bcScale = 1.0/h;
        bc = EssentialBC(left, bcScale*vPsi*(psi - control), bcQuad);
      }
    else
      {
        eqn = eqn + Integral(left, (1.0/h)*gammaBC*vPsi*(psi - control), bcQuad);
      }


    /* Create a discrete space, and discretize the function 0.0 on it */
    DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
    Expr psi0 = new DiscreteFunction(discSpace, 0.0, "psi0");

    /* create the forward problem */
    NonlinearProblem prob(mesh, eqn, bc, vPsi, psi, psi0, vecType);

    /* Write the sensitivity problems by hand. This will be unnecessary once
     * parametric differentiation is online. */
    Array<LinearProblem> sensProb(param.size());
    for (unsigned int i=0; i<param.size(); i++)
      { 
        Expr w = sin((i+1)*pi*y);
        Expr sensEqn = Integral(interior, (grad*vPsi)*(D*(grad*psi) - u*psi), quad)
          + Integral(interior, vPsi * psi * (alpha + 3.0*beta * psi0*psi0), quad)
          + Integral(right, vPsi*ux*psi0, bcQuad); 
        Expr sensBC;
        if (!robinBC)
          {
            sensBC = EssentialBC(left, bcScale*vPsi*(psi - w), bcQuad);
          }
        else
          {
            eqn = eqn + Integral(left, (1.0/h)*gammaBC*vPsi*(psi - w), bcQuad);
          }
        sensProb[i] = LinearProblem(mesh, sensEqn, sensBC, vPsi, psi, vecType);
      }


    /* Define the target expression */
    target_ = new DiscreteFunction(discSpace, 1.0);
    
    Expr objective = Integral(interior, 0.5*pow(psi - target_, 2.0), quad);
    Functional obj(mesh, objective, vecType);

    if (myRank==0) cout << "initializing NLP..." << endl;    

    initialize(param, psi, psi0, prob, sensProb, obj);

    

    if (myRank==0) cout << "done initializing NLP..." << endl; 
  }

  void CDRControlModel::setTargetVector(const Expr& expr)
  {
    const Vector<double>& vec 
      = DiscreteFunction::discFunc(expr)->getVector().copy();
    DiscreteFunction::discFunc(target_)->setVector(vec);
  }

}
#endif

