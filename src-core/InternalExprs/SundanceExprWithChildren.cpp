/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExprWithChildren.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


ExprWithChildren::ExprWithChildren(const Array<RefCountPtr<ScalarExpr> >& children)
	: EvaluatableExpr(), 
    children_(children),
    childDerivSetIndices_(children.size())
{}

int ExprWithChildren::setupEval(const RegionQuadCombo& region,
                                const EvaluatorFactory* factory) const
{

  /* If we've been here already with this region, or if we've
   * been here in another region having the same set of required
   * derivatives, we're done. */
     
  if (checkForKnownRegion(region))
    {
      return getDerivSetIndex(region);
    }

  /* If we have reached this line, the region and its required derivatives
   * are new. Create a new entry in our tables of deriv sets. 
   * This step also creates
   * a sparsity pattern for the current deriv set. */
  int derivSetIndex = registerRegion(region, currentDerivSuperset(), 
                                       factory);

  /* set up the children */
  childDerivSetIndices_.append(Array<int>());
  for (int i=0; i<children_.size(); i++)
    {
      int childDSI = evaluatableChild(i)->setupEval(region, factory);
      childDerivSetIndices_[i].append(childDSI);
    }

  return derivSetIndex;
  
}

bool ExprWithChildren::isConstant() const
{
  for (int i=0; i<children_.size(); i++) 
    {
      if (!children_[i]->isConstant()) return false;
    }
  return true;
}

void ExprWithChildren::accumulateUnkSet(Set<int>& unkIDs) const
{
  for (int i=0; i<children_.size(); i++) 
    {
      children_[i]->accumulateUnkSet(unkIDs);
    }
}

void ExprWithChildren::accumulateTestSet(Set<int>& testIDs) const
{
  for (int i=0; i<children_.size(); i++) 
    {
      children_[i]->accumulateTestSet(testIDs);
    }
}

const EvaluatableExpr* ExprWithChildren::evaluatableChild(int i) const
{
  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(children_[i].get());

  TEST_FOR_EXCEPTION(e==0, InternalError, 
                     "ExprWithChildren: cast of child [" 
                     << children_[i]->toString()
                     << " to evaluatable expr failed");

  return e;
}

void ExprWithChildren::getRoughDependencies(Set<Deriv>& funcs) const
{
  for (int i=0; i<children_.size(); i++)
    {
      evaluatableChild(i)->getRoughDependencies(funcs);
    }
}

void ExprWithChildren::resetReferenceCount() const
{
  refCount()++;

  /* Only update the reference count of the operands if this is the
   * first time this expression has been referenced. If this expression
   * is referenced multiple times, we still only need to evaluate the
   * operands once */
  if (refCount() == 1) 
    {
      for (int i=0; i<children_.size(); i++)
        {
          evaluatableChild(i)->resetReferenceCount();
        }
    }
}

void ExprWithChildren::flushResultCache() const
{
  resultCacheIsValid() = false;
  for (int i=0; i<children_.size(); i++)
    {
      evaluatableChild(i)->flushResultCache();
    }
}

bool ExprWithChildren::hasWorkspace() const
{
  TEST_FOR_EXCEPTION(refCount() <= 0, InternalError,
                     "bad ref count in ExprWithChildren::hasWorkspace()");
  refCount()--;
  if (refCount()==0) 
    {
      for (int i=0; i<children_.size(); i++)
        {
          evaluatableChild(i)->hasWorkspace();
        }
      return true;
    }
  return false;
}

void ExprWithChildren::resetDerivSuperset() const 
{
  currentDerivSuperset() = DerivSet();
  
  for (int i=0; i<children_.size(); i++)
    {
      evaluatableChild(i)->resetDerivSuperset();
    }
}

void ExprWithChildren::findDerivSuperset(const DerivSet& derivs) const 
{
  Tabs tabs;
  
  if (verbosity() > 1)
    {
      cerr << tabs << "finding deriv superset for expr " << toString()
           << endl;
      {
        Tabs tabs1;
        cerr << tabs1 << "required derivs are " << derivs.toString() << endl;
        cerr << tabs1 << "current derivs superset is " 
             << currentDerivSuperset().toString() << endl;
      }
    }
  
  currentDerivSuperset().merge(derivs);

  if (verbosity() > 1)
    {
      Tabs tabs1;
      cerr << tabs1 << "merged superset is " 
           << currentDerivSuperset().toString() << endl;
    }

  Array<DerivSet> childDerivs = derivsRequiredFromOperands(derivs);
  for (int i=0; i<children_.size(); i++)
    {
      Tabs tabs1;
      
      if (verbosity() > 1)
        {
          cerr << tabs1 << "recursing to child " << i 
               << ": " << children_[i]->toString() << endl;
        } 
      if (verbosity() > 1)
        {
          cerr << tabs1 << "child " << i << " requires " 
               << childDerivs[i].toString() << endl;
        }
      //      childDerivs[i].merge(currentDerivSuperset());
      evaluatableChild(i)->findDerivSuperset(childDerivs[i]);
    }
}

Array<DerivSet> ExprWithChildren::derivsRequiredFromOperands(const DerivSet& d) const
{
  Tabs tabs;

  if (verbosity() > 1)
    {
      cerr << tabs << "ExprWithChildren::derivsRequiredFromOperands()" << endl;
    }
  
  Array<DerivSet> rtn(children_.size());
  for (int i=0; i<children_.size(); i++)
    {
      rtn[i] = d;
    }
  return rtn;
}

