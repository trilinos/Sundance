/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

Evaluator::Evaluator()
{;}

CoordExprEvaluator::CoordExprEvaluator(const CoordExpr* expr)
  : SubtypeEvaluator<CoordExpr>(expr)
{;}

ConstantEvaluator::ConstantEvaluator(const SpatiallyConstantExpr* expr)
  : SubtypeEvaluator<SpatiallyConstantExpr>(expr)
{;}

SymbolicFuncElementEvaluator
::SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr)
  : SubtypeEvaluator<SymbolicFuncElement>(expr)
{;}

DiscreteFuncElementEvaluator
::DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr)
  : SubtypeEvaluator<DiscreteFuncElement>(expr)
{;}

void CoordExprEvaluator::eval(const EvalManager& mgr,
                              RefCountPtr<EvalVectorArray>& results) const 
{
  Tabs tabs;

  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---CoordExprEvaluator---");
  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const RefCountPtr<SparsityPattern>& s = expr()->sparsity(derivSetIndex);


  if (verbosity() > 1)
    {
      cerr << tabs << "CoordExprEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      cerr << tabs << "sparsity = " << endl << *s << endl;
    }

  results = mgr.stack().popVectorArray(s.get());

  EvalVectorArray& r = *results;

  for (int i=0; i<r.length(); i++) 
    {
      if (s->isZero(i)) r[i]->setToZero();
      else if (s->isConstant(i)) r[i]->setToOne();
      else 
        {
          mgr.evalCoordExpr(expr(), r[i]);
        }
    }
}

void SymbolicFuncElementEvaluator::eval(const EvalManager& mgr,
                                    RefCountPtr<EvalVectorArray>& results) const 
{
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---SymbolicFuncElementEvaluator---");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const RefCountPtr<SparsityPattern>& s = expr()->sparsity(derivSetIndex);

  results = mgr.stack().popVectorArray(s.get());

  if (verbosity() > 1)
    {
      cerr << tabs << "SymbolicFuncElementEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      cerr << tabs << "deriv set index = " << derivSetIndex << endl;
      cerr << tabs << "sparsity = " << endl << *s << endl;
    }

  EvalVectorArray& r = *results;

  for (int i=0; i<r.length(); i++) 
    {
      if (s->isZero(i)) r[i]->setToZero();
      else
        {
          switch(s->deriv(i).order())
            {
            case 0:
              {
                const ZeroExpr* z 
                  = dynamic_cast<const ZeroExpr*>(expr()->evalPt());
                if (z != 0)
                  {
                    r[i]->setToZero();
                    break;
                  }
                const DiscreteFuncElement* df 
                  = dynamic_cast<const DiscreteFuncElement*>(expr()->evalPt());
                TEST_FOR_EXCEPTION(df==0,
                                   InternalError,
                                   "SymbolicFuncElementEvaluator::eval() detected an"
                                   " unknown function with an evaluation point "
                                   "that is not a discrete function. The "
                                   "function is " << expr()->toString() 
                                   << ". Its eval point is " 
                                   << expr()->evalPt()->toString());
                mgr.evalDiscreteFuncElement(df, MultiIndex(), r[i]);
              }
              break;
            case 1: 
              r[i]->setToOne();
              break;
            default:
              TEST_FOR_EXCEPTION(true,
                                 InternalError,
                                 "SymbolicFuncElementEvaluator::eval() detected a "
                                 "structurally nonzero derivative of higher "
                                 "than first order");
            }
        }
      if (verbosity() > 1)
        {
          cerr << tabs << i << " " << s->deriv(i) << " " << r[i]->getStringValue()
               << endl;
        }
    }
}

void DiscreteFuncElementEvaluator::eval(const EvalManager& mgr,
                                     RefCountPtr<EvalVectorArray>& results) const 
{
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---DiscreteFuncElementEvaluator---");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const RefCountPtr<SparsityPattern>& s = expr()->sparsity(derivSetIndex);

  results = mgr.stack().popVectorArray(s.get());

  EvalVectorArray& r = *results;

  if (verbosity() > 1)
    {
      cerr << tabs << "DiscreteFuncElementEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      cerr << tabs << "sparsity = " << endl << *s << endl;
    }

  for (int i=0; i<r.length(); i++) 
    {
      if (s->isZero(i)) r[i]->setToZero();
      else
        {
          MultiIndex mi;
          if (s->isFirstOrderSpatialDeriv(i))
            {
              mi[s->spatialDerivDir(i)] = 1;
            }
          mgr.evalDiscreteFuncElement(expr(), mi, r[i]);
        }
      if (verbosity() > 1)
        {
          cerr << tabs << i << " " << s->deriv(i) << " " << r[i]->getStringValue()
               << endl;
        }
    }
}

void ConstantEvaluator::eval(const EvalManager& mgr,
                             RefCountPtr<EvalVectorArray>& results) const 
{
  Tabs tabs;
  SUNDANCE_OUT(verbosity() > VerbLow, tabs << "---ConstantEvaluator---");

  int derivSetIndex = expr()->getDerivSetIndex(mgr.getRegion());
  const RefCountPtr<SparsityPattern>& s = expr()->sparsity(derivSetIndex);

  results = mgr.stack().popVectorArray(s.get());

  EvalVectorArray& r = *results;

  if (verbosity() > 1)
    {
      cerr << tabs << "ConstantEvaluator::eval: expr=" << expr()->toString() << endl;
      cerr << tabs << "sparsity = " << endl << *s << endl;
    }

  for (int i=0; i<r.length(); i++) 
    {
      TEST_FOR_EXCEPTION(!(s->isConstant(i) || s->isZero(i)), InternalError,
                         "non-constant sparsity table "
                         "entry detected in ConstantEvaluator::eval()");
      if (s->isZero(i)) r[i]->setToZero();
      else r[i]->setToConstantValue(expr()->value());
      if (verbosity() > 1)
        {
          cerr << tabs << i << " " << s->deriv(i) << " " << r[i]->getStringValue()
               << endl;
        }
    }
}
