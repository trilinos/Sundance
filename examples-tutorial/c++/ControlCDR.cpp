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
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"

#include "SundanceNLPModelEvaluator.hpp"
#include "Sundance.hpp"



/** 
 * Optimal target matching for steady-state convection-diffusion equation with a cubic
 * reaction term.
 */

/* These macros define predicate classes to be used in identifying
 * boundary regions */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;});
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;});
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;});


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
    Vector<double> getInitialParameters() const 
    {
      Vector<double> rtn = paramSpace().createMember();
      for (int i=0; i<paramSpace().dim(); i++) rtn[i] = 0.25;
      return rtn;
    }

    const string& outfile() const {return outfile_;}

    const Mesh& mesh() const {return mesh_;}

  private:
    string outfile_;
    Mesh mesh_;
  };
}


int main(int argc, void** argv)
{
  using MoochoPack::MoochoSolver;
	using NLPInterfacePack::NLPFirstOrderThyraModelEvaluator;

  try
		{
      bool doSim=false;
      Sundance::setOption("doSim", "doOpt", doSim, "whether to do simulation or optimization");
      string solverSpec = "belos";
      Sundance::setOption("solver", solverSpec, "specify solver: options are amesos, aztec, belos");

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

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      RefCountPtr<Thyra::CDRControlModel> ssModel 
        = rcp(new Thyra::CDRControlModel(params, vecType));
      RefCountPtr<Thyra::ModelEvaluator<double> > model = ssModel;

      RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >lowsFactory;

      if (solverSpec=="aztec")
        {
          lowsFactory = rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
          lowsFactory->setPreconditionerFactory(rcp(new Thyra::IfpackPreconditionerFactory()),"");
        }
      else if (solverSpec=="belos")
        {
          lowsFactory = rcp(new Thyra::BelosLinearOpWithSolveFactory<double>());
          lowsFactory->setPreconditionerFactory(rcp(new Thyra::IfpackPreconditionerFactory()),"");
        }
      else
        {
          lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());
        }
      RefCountPtr<Thyra::ModelEvaluator<double> > modelWithLOWS
        = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model,
                                                                       lowsFactory));

      int flag=0;
      if (doSim)
        {
          flag = -1;
        }

      NLPFirstOrderThyraModelEvaluator nlp(modelWithLOWS, flag, flag);

      // Create the solver object
      MoochoSolver  solver;
      
      // Set the NLP
      solver.set_nlp( Teuchos::rcp(&nlp,false) );

      if (myRank==0) cout << "starting solve" << endl;
      // Solve the NLP
      const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();

      Expr uFinal = ssModel->stateVariable();
      FieldWriter writer = new VTKWriter(ssModel->outfile());
      writer.addMesh(ssModel->mesh());
      writer.addField("Optimal Psi", new ExprFieldWrapper(uFinal[0]));
      writer.write();
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
    : SundanceNLPModelEvaluator(vecType), outfile_(), mesh_()
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
    
    CellFilter left = edges.subset(new LeftPointTest());
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
    
    /* Create parameter expressions for alpha and beta. These will be
     * modified during the continuation loop, and thus must be
     * parameter exprs rather than constants. */
    Expr alpha = new Parameter(alpha0);
    Expr beta = new Parameter(beta0);

    
    /* Define the boundary field in terms of the design parameters */
    Expr control = 0.0;
    Array<Expr> p(numControls);
    for (unsigned int i=1; i<=p.size(); i++) 
      {
        p[i-1] = new Parameter(1.0);
        control = control + pi*pi*i*i*p[i-1] * sin(i*pi*y);
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
        
    /* Define the Dirichlet BC using the control field */
    Expr bc = EssentialBC(left, vPsi*(psi - control), bcQuad);


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
        Expr w = pi*pi*(i+1)*(i+1)*sin((i+1)*pi*y);
        Expr sensEqn = Integral(interior, (grad*vPsi)*(D*(grad*psi) - u*psi), quad)
          + Integral(interior, vPsi * psi * (alpha + 3.0*beta * psi0*psi0), quad)
          + Integral(right, vPsi*ux*psi0, bcQuad); 
        Expr sensBC = EssentialBC(left, vPsi*(psi - w), bcQuad);
        sensProb[i] = LinearProblem(mesh, sensEqn, sensBC, vPsi, psi, vecType);
      }

    Expr psiStar = sin(2.0*pi*y);// + 0.1*cos(5.0*pi*y);
    
    Expr objective = Integral(interior, 0.5*pow(psi - psiStar, 2.0)*pow(x,6.0), quad);
    Functional obj(mesh, objective, vecType);

    if (myRank==0) cout << "initializing NLP..." << endl;    

    initialize(param, psi, psi0, prob, sensProb, obj);

    if (myRank==0) cout << "done initializing NLP..." << endl; 
  }

}
#endif

