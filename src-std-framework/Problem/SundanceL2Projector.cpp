/* @HEADER@ */
/* @HEADER@ */

#include "SundanceL2Projector.hpp"
#include "TSFAztecSolver.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceGaussianQuadrature.hpp"

#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"

using namespace SundanceStdFwk;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

L2Projector::L2Projector(const DiscreteSpace& space, 
                         const Expr& expr)
  : prob_(), solver_()
{
  TEST_FOR_EXCEPTION(space.basis().size() != expr.size(),
                     RuntimeError,
                     "mismatched vector structure between basis and expr");
  
  TEST_FOR_EXCEPTION(space.basis().size() == 0,
                     RuntimeError,
                     "Empty basis?");
  
  Expr v = new TestFunction(space.basis()[0]);
  Expr u = new UnknownFunction(space.basis()[0]);
  
  for (int i=1; i<space.basis().size(); i++)
    {
      v.append(new TestFunction(space.basis()[i]));
      u.append(new UnknownFunction(space.basis()[i]));
    }

  CellFilter interior = new MaximalCellFilter();

  Expr eqn = Integral(interior, v*(u-expr), new GaussianQuadrature(4));
  Expr bc;

  prob_ = LinearProblem(space.mesh(), eqn, bc, v, u, space.vecType());
  

  /* Create an Aztec solver for solving the linear subproblems */
  std::map<int,int> azOptions;
  std::map<int,double> azParams;
  
  azOptions[AZ_solver] = AZ_cg;
  azOptions[AZ_precond] = AZ_dom_decomp;
  azOptions[AZ_subdomain_solve] = AZ_icc;
  azOptions[AZ_graph_fill] = 1;
  azOptions[AZ_max_iter] = 1000;
  azParams[AZ_tol] = 1.0e-13;
  
  solver_ = new AztecSolver(azOptions,azParams);
}



