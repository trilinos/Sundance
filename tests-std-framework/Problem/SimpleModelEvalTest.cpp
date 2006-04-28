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


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})


namespace Thyra
{
  class SimpleSundanceModel : public SundanceNLPModelEvaluator
  {
  public:
    /** */
    SimpleSundanceModel(const VectorType<double>& vecType);

    /** */
    Array<double> exactSoln() const ;

    /** */
    Mesh mesh() const {return mesh_;}
    
  private:
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
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      RefCountPtr<Thyra::SimpleSundanceModel> ssModel 
        = rcp(new Thyra::SimpleSundanceModel(vecType));
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
      
      // Solve the NLP
      const MoochoSolver::ESolutionStatus	solution_status = solver.solve_nlp();

      Array<double> exact = ssModel->exactSoln();
      Array<double> soln = ssModel->parameters();
      cout << "exact solution: " << exact << endl;
      cout << "numerical solution: " << soln << endl;
      double error = 0.0;
      for (unsigned int i=0; i<exact.size(); i++) error += pow(exact[i]-soln[i], 2.0);
      error = sqrt(error);

      bool OK = solution_status == MoochoSolver::SOLVE_RETURN_SOLVED;
      Sundance::passFailTest("MOOCHO converged", OK, error, 1.0e-3);
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}


namespace Thyra
{
  SimpleSundanceModel::SimpleSundanceModel(const VectorType<double>& vecType)
    : SundanceNLPModelEvaluator(vecType),
      mesh_()
  {
    MeshType meshType = new BasicSimplicialMeshType();
    int np = MPIComm::world().getNProc();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10, meshType);
    mesh_ = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
    CellFilter points = new DimensionalCellFilter(0);
    CellFilter leftPoint = points.subset(new LeftPointTest());
    CellFilter rightPoint = points.subset(new RightPointTest());

    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    Expr u = new UnknownFunction(new Lagrange(2), "u");
    Expr v = new TestFunction(new Lagrange(2), "v");

    /* Create differential operator and coordinate function */
    Expr dx = new SundanceCore::Derivative(0);
    Expr x = new CoordExpr(0);

    const double pi = 4.0*atan(1.0);

    Expr source = 0.0;
    Array<Expr> a(4);

    for (unsigned int i=1; i<=a.size(); i++) 
      {
        a[i-1] = new Parameter(1.0);
        source = source + pi*pi*i*i* a[i-1] * sin(i*pi*x);
      }
    Expr alpha = new ListExpr(a);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);

      
    /* Define the weak form */
    Expr eqn = Integral(interior, (dx*v)*(dx*u) + v*source, quad);
    /* Define the Dirichlet BC */
    Expr bc = EssentialBC(leftPoint, v*u, quad)
      + EssentialBC(rightPoint, v*u, quad);

    /* Create a discrete space, and discretize the function 1.0 on it */
    DiscreteSpace discSpace(mesh_, new Lagrange(2), vecType);
    Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

    /* create the forward problem */
    NonlinearProblem prob(mesh_, eqn, bc, v, u, u0, vecType);
      
    /* Write the sensitivity problems by hand. This will be unnecessary once
     * parametric differentiation is online. */
    Array<LinearProblem> sensProb(alpha.size());
    for (unsigned int i=0; i<alpha.size(); i++)
      { 
        double w = (i+1)*(i+1)*pi*pi;
        Expr sensEqn = Integral(interior, 
                                (dx*v)*(dx*u) + w*v*sin((i+1)*pi*x), quad);
        Expr sensBC = EssentialBC(leftPoint, v*u, quad)
          + EssentialBC(rightPoint, v*u, quad);
        sensProb[i] = LinearProblem(mesh_, sensEqn, sensBC, v, u, vecType);
      }

    Expr sourceBasis = List(sin(pi*x), sin(2.0*pi*x), sin(3.0*pi*x), sin(4.0*pi*x));
    Expr uStar = 0.0;
    for (unsigned int i=0; i<exactSoln().size(); i++)
      {
        uStar = uStar - sourceBasis[i] * exactSoln()[i];
      }
    
    Expr objective = Integral(interior, 0.5*pow(u - uStar, 2.0), quad);
    Functional obj(mesh_, objective, vecType);

    initialize(alpha, u, u0, prob, sensProb, obj);
  }

  Array<double> SimpleSundanceModel::exactSoln() const
  {
    return tuple(1.0, 0.3, 0.2, 0.1);
  }
}
#endif

