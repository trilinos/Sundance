/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_L2PROJECTOR_H
#define SUNDANCE_L2PROJECTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceLinearProblem.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

  /**
   * L2Projector projects an expression onto a DiscreteSpace. 
   */
  class L2Projector
  {
  public:
    /** */
    L2Projector(const DiscreteSpace& space, 
                const Expr& expr);

    /** */
    Expr project() const {return prob_.solve(solver_);}

  private:
    LinearProblem prob_;

    LinearSolver<double> solver_;
    
  };
}


#endif
