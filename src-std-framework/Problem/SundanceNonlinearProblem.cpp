/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNonlinearProblem.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


NonlinearProblem::NonlinearProblem() 
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(),
    discreteU0_(0)
{}


NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
                                   const Expr& eqn, 
                                   const Expr& bc,
                                   const Expr& test, 
                                   const Expr& unk, 
                                   const Expr& u0, 
                                   const VectorType<double>& vecType)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    discreteU0_(0)
    
{
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, test, unk, u0, 
                          rcp(new BruteForceEvaluatorFactory())));

  assembler_ = rcp(new Assembler(mesh, eqnSet, vecType));

  discreteU0_ = dynamic_cast<DiscreteFunction*>(u0_.ptr().get());

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem ctor");

  VectorSpace<double> domain = assembler_->solutionSpace()->vecSpace();
  VectorSpace<double> range = assembler_->rowSpace()->vecSpace();

  setDomainAndRange(domain, range);
}

TSFExtended::Vector<double> NonlinearProblem::getInitialGuess() const 
{
  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::getInitialGuess()");

  return discreteU0_->vector();
}



LinearOperator<double> NonlinearProblem::
computeJacobianAndFunction(Vector<double>& functionValue) const
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEST_FOR_EXCEPTION(discreteU0_==0, RuntimeError,
                     "null discrete function pointer in "
                     "NonlinearProblem::jacobian()");

  TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, RuntimeError,
                     "null evaluation point in "
                     "NonlinearProblem::jacobian()");

  discreteU0_->setVector(currentEvalPt());

  assembler_->assemble(J_, functionValue);

  return J_;
}


