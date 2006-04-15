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
#include "TSF_config.h"

#ifdef HAVE_ENABLED_MOOCHO

#include "MoochoPack_MoochoSolver.hpp"
#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"

#include "Sundance.hpp"
#include "SundanceModelEvaluatorBase.hpp"

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;});
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;});


namespace Thyra
{
  class SimpleSundanceModel : public SundanceModelEvaluator
  {
  public:
    SimpleSundanceModel(const VectorType<double>& vecType)
      : SundanceModelEvaluator(vecType),
        paramSpace_(rcp(new DefaultSerialVectorSpace<double>(4))),
        alpha_(4),
        u0_(),
        prob_(),
        sensProb_(4),
        obj_(),
        objEval_()
    {
      MeshType meshType = new BasicSimplicialMeshType();
      int np = MPIComm::world().getNProc();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);
      Mesh mesh = mesher.getMesh();

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
      
      for (unsigned int i=1; i<=alpha_.size(); i++) 
        {
          alpha_[i-1] = new Parameter(1.0);
          source = source + pi*pi*i*i* alpha_[i-1] * sin(i*pi*x);
        }

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      
      /* Define the weak form */
      Expr eqn = Integral(interior, (dx*v)*(dx*u) + v*source, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*u, quad)
        + EssentialBC(rightPoint, v*u, quad);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      u0_ = new DiscreteFunction(discSpace, 1.0, "u0");

      /* create the forward problem */
      prob_  = rcp(new NonlinearProblem(mesh, eqn, bc, v, u, u0_, vecType));

      /* Write the sensitivity problems by hand. This will be unnecessary once
       * parametric differentiation is online. */
      for (unsigned int i=0; i<alpha_.size(); i++)
        { 
          double w = (i+1)*(i+1)*pi*pi;
          Expr sensEqn = Integral(interior, (dx*v)*(dx*u) + w*v*sin((i+1)*pi*x), quad);
          Expr sensBC = EssentialBC(leftPoint, v*u, quad)
            + EssentialBC(rightPoint, v*u, quad);
          sensProb_[i] = LinearProblem(mesh, sensEqn, sensBC, v, u, vecType);
        }

      Expr uStar = sin(pi*x) + 0.3*sin(2.0*pi*x);
      Expr objective = Integral(interior, 0.5*pow(u - uStar, 2.0), quad);
      obj_ = Functional(mesh, objective, vecType);
      objEval_ = obj_.evaluator(u, u0_);
    }

    /** */
    void internalEvalModel(const Vector<double>& stateVec,
                           const Vector<double>& params,
                           Vector<double>& resid,
                           double& objFuncVal,
                           LinearOperator<double>& df_dx,
                           Array<Vector<double> >& df_dp,
                           Vector<double>& dg_dp_T,
                           Vector<double>& dg_dx_T) const 
    {
      TEST_FOR_EXCEPTION(params.ptr().get()==0, RuntimeError,
                         "null params vector!");
      TEST_FOR_EXCEPTION( (int) alpha_.size() != params.space().dim(),
                         RuntimeError,
                         "Mismatch between input parameter vector space and "
                         "symbolic parameter object size");

      TEST_FOR_EXCEPTION(stateVec.ptr().get()==0, RuntimeError,
                         "null state vector!");


      const DefaultSerialVector<double>* dsv_params 
        = dynamic_cast<const DefaultSerialVector<double>*>(params.ptr().get());
      //      DefaultSerialVector<double>* dsv_dg_dp_T 
      // = dynamic_cast<DefaultSerialVector<double>*>(dg_dp_T.ptr().get());

      TEST_FOR_EXCEPTION(dsv_params == 0, RuntimeError, "params vector is not a default serial vector");
      //      TEST_FOR_EXCEPTION(dsv_dg_dp_T == 0, RuntimeError, "dg_dp_T vector is not a default serial vector");

      /* set the symbolic parameter values to the input parameters */
      cout << "parameters are: " << endl;
      for (unsigned int i=0; i<alpha_.size(); i++)
        {

          alpha_[i].setParameterValue(dsv_params->getPtr()[i]);
          cout << i << " " << alpha_[i] << endl;
        }
      

      cout << "printing state vector" << endl;
      /* set the symbolic state function value to the input state */
      //      DiscreteFunction::discFunc(u0_)->setVector(stateVec);

      cout << "state vector is " << stateVec.description() << endl;
      prob_->setEvalPt(stateVec);

      /* compute the Jacobian and residual of the constraints */
      if (df_dx.ptr().get()!=0 && resid.ptr().get() != 0)
        {
          cout << "computing Jacobian and residual..." << endl;
          prob_->computeJacobianAndFunction(df_dx, resid);          
          cout << "residual vector is " << resid.description() << endl;
          cout << "J is " << endl;
          df_dx.print(cout);
        }
      else if (resid.ptr().get() != 0 && df_dx.ptr().get()==0)
        {
          cout << "computing residual w/o Jacobian..." << endl;
          prob_->computeFunctionValue(resid);
          cout << "residual vector is " << resid.description() << endl;
        }
      else if (df_dx.ptr().get()!=0 && resid.ptr().get() == 0)
        {
          cout << "computing Jacobian w/o residual..." << endl;
          Vector<double> dummy = constraintSpace().createMember();
          prob_->computeJacobianAndFunction(df_dx, dummy);          
          cout << "J is " << endl;
          df_dx.print(cout);
        }
      else
        {
          cout << "requested neither the residual nor the Jacobian" << endl;
        }
        

      /* compute the derivatives of the constraint residual wrt the parameters */
      cout << "sensitivities are: " << endl;
      for (unsigned int i=0; i<alpha_.size(); i++)
        {
          df_dp[i] = sensProb_[i].getRHS();
          cout << i << " " << df_dp[i].description() << endl;
        }

      /* evaluate the objective function and its gradient wrt the states */ 
      if (dg_dx_T.ptr().get() != 0)
        {
          cout << "computing objective function and gradient" << endl;
          Expr df_dx_Expr = objEval_.evalGradient(objFuncVal);
          dg_dx_T.acceptCopyOf(DiscreteFunction::discFunc(df_dx_Expr)->getVector());
          cout << "value = " << objFuncVal << endl;
          cout << "gradient = " << dg_dx_T.description() << endl;
        }
      else
        {
          cout << "computing objective function" << endl;
          objFuncVal = objEval_.evaluate();
          cout << "value = " << objFuncVal << endl;
        }
      


      /* evaluate the derivatives of the objective function wrt the parameters */
      /* set the symbolic parameter values to the input parameters */
      if (dg_dp_T.ptr().get() != 0)
        {
          dg_dp_T.setToConstant(0.0);
        }
    }

    /** */
    VectorSpace<double> paramSpace() const 
    {
      return paramSpace_;
    }        
           
    /** */
    VectorSpace<double> stateSpace() const 
    {
      VectorSpace<double> rtn = createW().domain();
      cout << "state space is " << rtn.description() << endl;
      return rtn;
    }
           
    /** */
    VectorSpace<double> constraintSpace() const 
    {
      VectorSpace<double> rtn = createW().range();
      cout << "constraint space is " << rtn.description() << endl;
      return rtn;
    }
           
    /** */
    LinearOperator<double> createW() const 
    {
      static LinearOperator<double> J = prob_->allocateJacobian();
      TEST_FOR_EXCEPTION(J.ptr().get()==0, RuntimeError,
                         "null Jacobian");
      return J;
    }


    /** */
    Vector<double> getInitialState() const 
    {
      return stateSpace().createMember();
    }

    
    /** */
    Vector<double> getInitialParameters() const 
    {
      Vector<double> rtn = paramSpace().createMember();
      rtn.setToConstant(3.14159);
      return rtn;
    }

    
      

  private:
    VectorSpace<double> paramSpace_;

    mutable Array<Expr> alpha_;

    mutable Expr u0_;

    RefCountPtr<NonlinearProblem> prob_;

    Array<LinearProblem> sensProb_;

    Functional obj_;

    FunctionalEvaluator objEval_;

    

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
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      RefCountPtr<Thyra::ModelEvaluator<double> > model 
        = rcp(new Thyra::SimpleSundanceModel(vecType));

      RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
        lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());

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
      //#else
      cout << "moocho not present: test INACTIVE" << endl;

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
#else

#error blah

#include <iostream>

int main()
{
  std::cout << "moocho not present: test INACTIVE" << endl;
}

#endif
