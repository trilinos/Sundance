/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellJacobianBatch.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

static Time& functionalEvalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("functional evaluation"); 
  return *rtn;
}

static Time& functionalEvalSetupTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("functional evaluation setup"); 
  return *rtn;
}

namespace SundanceStdFwk
{
  double evaluateIntegral(const Mesh& mesh, const Expr& expr)
  {
    FunctionalEvaluator eval(mesh, expr);
    return eval.evaluate();
  }
}

FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral)
  : assembler_()
{
  RefCountPtr<EquationSet> eqnSet = rcp(new EquationSet(integral));

  assembler_ = rcp(new Assembler(mesh, eqnSet));
}


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral
                                         const Expr& bcs,
                                         const Expr& fields,
                                         const Expr& fieldValues)
  : assembler_()
{
  Expr f = fields.flatten();
  Expr f0 = fieldValues.flatten();
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0));

  assembler_ = rcp(new Assembler(mesh, eqnSet));
}

double FunctionalEvaluator::evaluate() const
{
  double value;
  assembler_->evaluate(value);
  return value;
}
