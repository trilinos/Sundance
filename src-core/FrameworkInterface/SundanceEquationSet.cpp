/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEquationSet.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

 


using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;

EquationSet::EquationSet(const Expr& eqns, 
                         const Expr& bcs,
                         const Expr& tests, 
                         const Expr& unks,
                         const Expr& unkLinearizationPts)
  : regions_(),
    testsOnRegions_(),
    unksOnRegions_(),
    testUnkPairsOnRegions_(),
    bcTestUnkPairsOnRegions_(),
    bcTestsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(2),
    bcRegionQuadComboNonzeroDerivs_(2),
    rqcToContext_(2),
    bcRqcToContext_(2),
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

  
  /* map each test and unknown function's ID numbers to its
   * position in the input function lists */
  for (int i=0; i<tests.size(); i++)
    {
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(tests[i].ptr().get());
      TEST_FOR_EXCEPTION(t==0, RuntimeError, 
                         "EquationSet ctor input test function "
                         << tests[i] 
                         << " does not appear to be a test function");
      int fid = t->funcID();
      testIDToReducedIDMap_.put(fid, i);
    }
  for (int i=0; i<unks.size(); i++)
    {
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(unks[i].ptr().get());
      TEST_FOR_EXCEPTION(u==0, RuntimeError, 
                         "EquationSet ctor input unk function "
                         << unks[i] 
                         << " does not appear to be a unk function");
      int fid = u->funcID();
      unkIDToReducedIDMap_.put(fid, i);
    }

  

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

  Array<int> contextID = tuple(EvalContext::nextID(),
                               EvalContext::nextID());

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
          for (int order=1; order<=2; order++)
            {
              EvalContext context(rqc, order, contextID[order-1]);
              DerivSet nonzeros 
                = SymbPreprocessor::setupExpr(term, tests,
                                              unks, 
                                              unkLinearizationPts,
                                              context);
              if (order==2) addToTestUnkPairs(rqc.domain(), nonzeros, false);
              rqcToContext_[order-1].put(rqc, context);
              regionQuadComboNonzeroDerivs_[order-1].put(rqc, nonzeros);
            }
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
          cerr << "working on region " << reg << endl;
          cerr << "instance ID= " << reg.ptr()->id() << endl;
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

          cerr << "BC regions are " << endl;
          for (Map<OrderedHandle<CellFilterStub>, Set<int> >::const_iterator 
                 iter=bcTestsOnRegions_.begin(); iter != bcTestsOnRegions_.end();
               iter ++)
            {
              cerr << "region=" << iter->first << endl;
              cerr << "id = " << iter->first.ptr()->id() << endl;
            }
          cerr << "checking map integrity " << endl;
          TEST_FOR_EXCEPTION(!bcTestsOnRegions_.containsKey(reg),
                             InternalError,
                             "Bug: region " << reg << " should appear in "
                             "BC list" << bcTestsOnRegions_);
          for (int t=0; t<bcSum->numTerms(d); t++)
            {
              RegionQuadCombo rqc(reg.ptr(), bcSum->quad(d,t));
              Expr term = bcSum->expr(d,t);
              rqcBCSet.put(rqc);
              bcRegionQuadComboExprs_.put(rqc, bcSum->expr(d,t)); 
              for (int order=1; order<=2; order++)
                {
                  EvalContext context(rqc, order, contextID[order-1]);
                  DerivSet nonzeros 
                    = SymbPreprocessor::setupExpr(term, tests,
                                                  unks, 
                                                  unkLinearizationPts,
                                                  context);
                  if (order==2) addToTestUnkPairs(rqc.domain(), nonzeros, true);
                  bcRqcToContext_[order-1].put(rqc, context);
                  bcRegionQuadComboNonzeroDerivs_[order-1].put(rqc, nonzeros);
                }
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

  for (int r=0; r<numRegions(); r++)
    {
      cerr << "region " << regions_[r] << endl
           << "isBCRegion=" << isBCRegion(r)
           << endl;
    }

}



void EquationSet
::addToTestUnkPairs(const OrderedHandle<CellFilterStub>& domain,
                    const DerivSet& nonzeros, 
                    bool isBC)
{
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "finding test-unk pairs "
               "for domain " << domain);
  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "isBC=" << isBC);
  
  RefCountPtr<Set<OrderedPair<int, int> > > funcPairs;
  Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > >* theMap;

  if (isBC) 
    {
      theMap = &(bcTestUnkPairsOnRegions_);
    }
  else 
    {
      theMap = &(testUnkPairsOnRegions_);
    } 
  SUNDANCE_OUT(verbosity() > VerbHigh, tab << "map ptr=" << theMap);
  if (theMap->containsKey(domain))
    {
      funcPairs = theMap->get(domain);
    }
  else
    {
      funcPairs = rcp(new Set<OrderedPair<int, int> >());
      theMap->put(domain, funcPairs);
    }

  for (DerivSet::const_iterator i=nonzeros.begin(); i!=nonzeros.end(); i++)
    {
      const MultipleDeriv& md = *i;
      if (md.order() != 2) continue;
      
      int testID = -1;
      int unkID = -1;

      MultipleDeriv::const_iterator j;
      for (j=md.begin(); j != md.end(); j++)
        {
          const Deriv& d = *j;
          const FunctionalDeriv* fd = d.funcDeriv();
          TEST_FOR_EXCEPTION(fd==0, InternalError, "non-functional deriv "
                             << d << " detected in EquationSet::"
                             "addToTestUnkPairs()");
          const UnknownFuncElement* u 
            = dynamic_cast<const UnknownFuncElement*>(fd->func());
          const TestFuncElement* t 
            = dynamic_cast<const TestFuncElement*>(fd->func());

          if (u != 0)
            {
              unkID = fd->funcID();
              continue;
            }
          if (t != 0)
            {
              testID = fd->funcID();
              continue;
            }
        }
      TEST_FOR_EXCEPTION(unkID==-1 || testID==-1, InternalError,
                         "multiple derivative " << md << " does not "
                         "appear to have exactly one test and unknown each");
      funcPairs->put(OrderedPair<int, int>(testID, unkID));
    }

  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "found " << *funcPairs);
  
}

const RefCountPtr<Set<OrderedPair<int, int> > >& EquationSet::
bcTestUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
{
  TEST_FOR_EXCEPTION(!bcTestUnkPairsOnRegions_.containsKey(domain),
                     InternalError,
                     "equation set does not have a test-unk pair list for "
                     "bc region " << domain);
  const RefCountPtr<Set<OrderedPair<int, int> > >& rtn 
    = bcTestUnkPairsOnRegions_.get(domain);

  TEST_FOR_EXCEPTION(rtn.get()==0, InternalError, 
                     "null test-unk pair list for BC region " << domain);
  return rtn;
}

bool EquationSet::isBCRegion(int d) const
{
  return bcTestsOnRegions_.containsKey(regions_[d]);
}



