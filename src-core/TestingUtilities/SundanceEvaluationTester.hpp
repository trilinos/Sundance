/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATIONTESTER_H
#define SUNDANCE_EVALUATIONTESTER_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTestEvalMediator.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceADCoord.hpp"
#include "SundanceADDerivative.hpp"
#include "TSFObjectWithVerbosity.hpp"


namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class EvaluationTester 
    : public TSFExtended::ObjectWithVerbosity<EvaluationTester>
  {
  public:
    /** */
    EvaluationTester(const Expr& e);

    /** */
    double evaluate(Array<double>& firstDerivs, 
                    Array<Array<double> >& secondDerivs) const ;
    
    double fdEvaluate(const double& step, const double& tol, 
                      const double& tol2,
                      bool& isOK);

    int numNonzeros() const {return sparsity_->numDerivs();}

    int numNodes() const {return ev_->countNodes();}

  private:
    
    Expr e_;
    RegionQuadCombo rqc_;
    EvalContext context_;
    EvalManager mgr_;
    RefCountPtr<AbstractEvalMediator> mediator_;
    mutable TestEvalMediator* tem_;
    const EvaluatableExpr* ev_;
    RefCountPtr<SparsitySuperset> sparsity_;
    Map<int, int> unkIDToDiscreteIDMap_;
  };

}



#endif
