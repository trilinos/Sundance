/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


SymbolicFuncElementEvaluator
::SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<SymbolicFuncElement>(expr, context),
    mi_(),
    spatialDerivs_(),
    ones_(),
    df_(dynamic_cast<const DiscreteFuncElement*>(expr->evalPt())),
    stringReps_()
{
  
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing symbolic func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  const ZeroExpr* z 
    = dynamic_cast<const ZeroExpr*>(expr->evalPt());
  
  TEST_FOR_EXCEPTION(z==0 && df_==0, InternalError,
                     "SymbolicFuncElementEvaluator ctor detected an "
                     "evaluation point=" << expr->toString()
                     << " that is neither zero nor a discrete "
                     "function.");

  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  
  int constantCounter = 0;
  int vectorCounter = 0;

  Set<MultiIndex> miSet;
  
  for (int i=0; i<sparsity()->numDerivs(); i++) 
    {
      if (sparsity()->isSpatialDeriv(i))
        {
          /* evaluate the spatial deriv applied to the evaluation point */
          TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "spatial derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");

          mi_.append(sparsity()->multiIndex(i));
          miSet.put(sparsity()->multiIndex(i));
          addVectorIndex(i, vectorCounter);
          spatialDerivs_.append(vectorCounter++);
          int dir = sparsity()->multiIndex(i).firstOrderDirection();
          string deriv = "D[" + df_->name() + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
      else
        {
          TEST_FOR_EXCEPTION(sparsity()->deriv(i).order() > 1,
                             InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "nonzero functional derivative of order greater "
                             "than one. All such derivs should have been "
                             "identified as zero by this point. The bad "
                             "derivative is " << sparsity()->deriv(i)
                             << ", and the bad sparsity table is "
                             << *sparsity());

          if (sparsity()->deriv(i).order()==0)
            {
              TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "zero-order derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");
              /* value of zeroth functional deriv is a discrete function */
              addVectorIndex(i, vectorCounter);
              spatialDerivs_.append(vectorCounter++);
              mi_.append(MultiIndex());
              miSet.put(MultiIndex());
              stringReps_.append(df_->name());
            }
          else
            {
              /* value of first functional deriv is one */
              addConstantIndex(i, constantCounter);
              ones_.append(constantCounter++);
            }
        }
    }

  if (df_ != 0)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "setting up evaluation for discrete eval pt");
      df_->setupEval(context);
      dfEval_ = dynamic_cast<const DiscreteFuncElementEvaluator*>(df_->evaluator(context).get());
    }
}




void SymbolicFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(symbolicFuncEvalTimer());
  Tabs tabs;
  
  if (verbosity() > VerbSilent)
    {
      cerr << tabs << "SymbolicFuncElementEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      if (verbosity() > VerbLow)
        {
          cerr << tabs << "sparsity = " << endl << *sparsity() << endl;
        }
    }

  constantResults.resize(ones_.size());
  vectorResults.resize(spatialDerivs_.size());

  /* Evaluate discrete functions if necessary */
  if (df_ != 0 && mi_.size() > 0)
    {
      for (unsigned int i=0; i<mi_.size(); i++)
        {
          vectorResults[i] = mgr.popVector();
          TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                             InternalError,
                             "invalid evaluation vector allocated in "
                             "SymbolicFuncElementEvaluator::internalEval()");
          vectorResults[i]->setString(stringReps_[i]);
        }
      mgr.evalDiscreteFuncElement(df_, mi_, vectorResults);
    }

  /* Set the known one entries to one */
  for (unsigned int i=0; i<ones_.size(); i++)
    {
      constantResults[ones_[i]] = 1.0;
    }
  

}


