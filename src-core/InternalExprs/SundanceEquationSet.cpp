/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEquationSet.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"

 

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;

EquationSet::EquationSet(const Expr& eqns, 
            const Expr& bcs,
            const Expr& tests, 
            const Expr& unks,
            const Expr& unkLinearizationPts)
{
  verbosity() = classVerbosity();

  /* begin with a sanity check to ensure that the input equation set 
   * exists and is integral form */
  SUNDANCE_OUT(verbosity() > VerbLow, 
               "checking existence of input eqn set...");
  const SumOfIntegrals* integralSum
    = dynamic_cast<const SumOfIntegrals*>(eqns.ptr().get());

  TEST_FOR_EXCEPTION(eqns.ptr().get()==0, RuntimeError,
                     "EquationSet ctor detected empty equation set input");

  TEST_FOR_EXCEPTION(integralSum==0, RuntimeError,
                     "EquationSet ctor detected an input equation set that "
                     "is not in integral form");
  SUNDANCE_OUT(verbosity() > VerbLow, 
               "...input eqn set is OK");

  /* determine whether or not this problem includes essential BCs */
  SUNDANCE_OUT(verbosity() > VerbLow, 
               "checking whether the eqn set includes essential BCs...");
  bool hasBCs = false;
  const SumOfBCs* bcSum = 0 ;
  if (bcs.ptr().get() != 0)
    {
      bcSum = dynamic_cast<const SumOfBCs*>(bcs.ptr().get());
      TEST_FOR_EXCEPTION(bcSum==0, RuntimeError,
                         "EquationSet ctor detected an input Essential "
                         "BC that is not an EssentialBC object. "
                         "Object is " << bcs);
      hasBCs = true;
    }
  if (hasBCs)
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "...problem has EssentialBCs");
    }
  else
    {
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   "...problem has no EssentialBCs");
    }

  /* Now compile a list of all domains appearing in either the eqns or
   * the BCs */
  for (int d=0; d<integralSum->numDomains(); d++)
    {
      OrderedHandle<CellFilterBase> dom = integralSum->domain(d);
      if (!domains_.contains(dom)) 
        {
          domains_.put(dom);
          testsOnDomains_.put(dom, integralSum->testsOnDomain(d));
          unksOnDomains_.put(dom, integralSum->unksOnDomain(d));
        }
      else
        {
          Set<int>& t = testsOnDomains_.get(dom);
          t.merge(integralSum->testsOnDomain(d));
            Set<int>& u = unksOnDomains_.get(dom);
          u.merge(integralSum->unksOnDomain(d));
        }
    }
  
  if (hasBCs)
    {
      /* functions found in the BCs both in the overall lists and 
       * also in the bc-specific lists */
      for (int d=0; d<bcSum->numDomains(); d++)
        {
          OrderedHandle<CellFilterBase> dom = bcSum->domain(d);
          if (!domains_.contains(dom)) 
            {
              domains_.put(dom);
              testsOnDomains_.put(dom, bcSum->testsOnDomain(d));
              unksOnDomains_.put(dom, bcSum->unksOnDomain(d));
              bcTestsOnDomains_.put(dom, bcSum->testsOnDomain(d));
              bcUnksOnDomains_.put(dom, bcSum->unksOnDomain(d));
            }
          else
            {
              if (!bcTestsOnDomains_.containsKey(dom))
                {
                  bcTestsOnDomains_.put(dom, bcSum->testsOnDomain(d));
                }
              if (!bcUnksOnDomains_.containsKey(dom))
                {
                  bcUnksOnDomains_.put(dom, bcSum->unksOnDomain(d));
                }
              Set<int>& t = testsOnDomains_.get(dom);
              t.merge(bcSum->testsOnDomain(d));
              Set<int>& u = unksOnDomains_.get(dom);
              u.merge(bcSum->unksOnDomain(d));
            }
        }
    }

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Tests appearing on each domain: " << endl << testsOnDomains_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each domain: " << endl << unksOnDomains_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Tests appearing on each BC domain: " 
               << endl << bcTestsOnDomains_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each BC domain: " 
               << endl << bcUnksOnDomains_);

}
