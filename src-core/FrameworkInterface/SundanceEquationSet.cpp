/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEquationSet.hpp"
#include "SundanceSymbPreprocessor.hpp"
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
                         const Expr& unkLinearizationPts,
                         const RefCountPtr<EvaluatorFactory>& evalFactory)
  : regions_(),
    testsOnRegions_(),
    unksOnRegions_(),
    bcTestsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    testFuncs_(tests),
    unkFuncs_(unks),
    unkLinearizationPts_(unkLinearizationPts),
    testIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    isNonlinear_(false)
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


  Set<OrderedHandle<CellFilterStub> > regionSet;
  Set<RegionQuadCombo> rqcSet;
  Set<RegionQuadCombo> rqcBCSet;

  /* Now compile a list of all regions appearing in either the eqns or
   * the BCs */

  /* Do the non-bc eqns first */
  for (int d=0; d<integralSum->numRegions(); d++)
    {
      OrderedHandle<CellFilterStub> reg = integralSum->region(d);
      if (!regionSet.contains(reg)) 
        {
          regionSet.put(reg);
          testsOnRegions_.put(reg, integralSum->testsOnRegion(d));
          unksOnRegions_.put(reg, integralSum->unksOnRegion(d));
        }
      else
        {
          Set<int>& t = testsOnRegions_.get(reg);
          t.merge(integralSum->testsOnRegion(d));
          Set<int>& u = unksOnRegions_.get(reg);
          u.merge(integralSum->unksOnRegion(d));
        }
      for (int t=0; t<integralSum->numTerms(d); t++)
        {
          RegionQuadCombo rqc(reg.ptr(), integralSum->quad(d,t));
          Expr term = integralSum->expr(d,t);
          rqcSet.put(rqc);
          regionQuadComboExprs_.put(rqc, term);
          DerivSet nonzeros = SymbPreprocessor::setupExpr(term, tests,
                                                          unks, 
                                                          unkLinearizationPts,
                                                          rqc,
                                                          evalFactory.get());
          regionQuadComboNonzeroDerivs_.put(rqc, nonzeros);
        }
    }
  
  /* now do the BCs */
  if (hasBCs)
    {
      /* functions found in the BCs both in the overall lists and 
       * also in the bc-specific lists */
      for (int d=0; d<bcSum->numRegions(); d++)
        {
          OrderedHandle<CellFilterStub> reg = bcSum->region(d);
          if (!regionSet.contains(reg)) 
            {
              regionSet.put(reg);
              testsOnRegions_.put(reg, bcSum->testsOnRegion(d));
              unksOnRegions_.put(reg, bcSum->unksOnRegion(d));
              bcTestsOnRegions_.put(reg, bcSum->testsOnRegion(d));
              bcUnksOnRegions_.put(reg, bcSum->unksOnRegion(d));
            }
          else
            {
              if (!bcTestsOnRegions_.containsKey(reg))
                {
                  bcTestsOnRegions_.put(reg, bcSum->testsOnRegion(d));
                }
              if (!bcUnksOnRegions_.containsKey(reg))
                {
                  bcUnksOnRegions_.put(reg, bcSum->unksOnRegion(d));
                }
              Set<int>& t = testsOnRegions_.get(reg);
              t.merge(bcSum->testsOnRegion(d));
              Set<int>& u = unksOnRegions_.get(reg);
              u.merge(bcSum->unksOnRegion(d));
            }
        
          for (int t=0; t<bcSum->numTerms(d); t++)
            {
              RegionQuadCombo rqc(reg.ptr(), bcSum->quad(d,t));
              Expr term = bcSum->expr(d,t);
              rqcBCSet.put(rqc);
              bcRegionQuadComboExprs_.put(rqc, bcSum->expr(d,t));
              DerivSet nonzeros 
                = SymbPreprocessor::setupExpr(term, tests,
                                              unks, 
                                              unkLinearizationPts,
                                              rqc,
                                              evalFactory.get());
              bcRegionQuadComboNonzeroDerivs_.put(rqc, nonzeros);
            }
        }
    }

  

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Tests appearing on each region: " << endl << testsOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each region: " << endl << unksOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Tests appearing on each BC region: " 
               << endl << bcTestsOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each BC region: " 
               << endl << bcUnksOnRegions_);


  /* convert sets to arrays */
  regions_ = regionSet.elements();
  regionQuadCombos_ = rqcSet.elements();
  bcRegionQuadCombos_ = rqcBCSet.elements();
  

}
