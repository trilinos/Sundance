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
                         const Expr& fields,
                         const Expr& fieldValues)
  : regions_(),
    varsOnRegions_(),
    unksOnRegions_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    varFuncs_(),
    unkFuncs_(),
    unkLinearizationPts_(),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true)
{
  Expr unks;
  Expr unkEvalPt;
  Expr vars;
  Expr varEvalPt;
  
  compTypes_.put(FunctionalOnly);

  rqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(FunctionalOnly,
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalOnly,
                                      Map<RegionQuadCombo, DerivSet>());


  init(eqns, bcs, vars, varEvalPt,
       unks, unkEvalPt,
       fields, fieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
                         const Expr& bcs,
                         const Expr& vars, 
                         const Expr& unks,
                         const Expr& unkLinearizationPts)
  : regions_(),
    varsOnRegions_(),
    unksOnRegions_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    varFuncs_(vars),
    unkFuncs_(unks),
    unkLinearizationPts_(unkLinearizationPts),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(false),
    isFunctionalCalculator_(false)
{
  Expr fixed;

  compTypes_.put(MatrixAndVector);
  compTypes_.put(VectorOnly);

  rqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
                                      Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(VectorOnly, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(VectorOnly, 
                                      Map<RegionQuadCombo, DerivSet>());


  init(eqns, bcs, vars, fixed,
       unks, unkLinearizationPts,
       fixed, fixed);
}


EquationSet::EquationSet(const Expr& eqns, 
                         const Expr& bcs,
                         const Expr& vars,
                         const Expr& varLinearizationPts, 
                         const Expr& unks,
                         const Expr& unkLinearizationPts,
                         const Expr& fixedFields,
                         const Expr& fixedFieldValues)
  : regions_(),
    varsOnRegions_(),
    unksOnRegions_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    varFuncs_(vars),
    unkFuncs_(unks),
    unkLinearizationPts_(unkLinearizationPts),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(false)
{

  compTypes_.put(MatrixAndVector);
  compTypes_.put(VectorOnly);

  rqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
                                      Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(VectorOnly, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(VectorOnly, 
                                      Map<RegionQuadCombo, DerivSet>());

  init(eqns, bcs, vars, varLinearizationPts, 
       unks, unkLinearizationPts,
       fixedFields, fixedFieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
                         const Expr& bcs,
                         const Expr& vars,
                         const Expr& varLinearizationPts, 
                         const Expr& fixedFields,
                         const Expr& fixedFieldValues)
  : regions_(),
    varsOnRegions_(),
    unksOnRegions_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    varFuncs_(vars),
    unkFuncs_(),
    unkLinearizationPts_(),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true)
{
  compTypes_.put(FunctionalOnly);
  compTypes_.put(FunctionalAndGradient);

  rqcToContext_.put(FunctionalAndGradient, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalAndGradient, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(FunctionalAndGradient, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalAndGradient, 
                                      Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(FunctionalOnly, 
                                    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalOnly, 
                                      Map<RegionQuadCombo, DerivSet>());

  init(eqns, bcs, vars, varLinearizationPts, 
       unkFuncs_, unkLinearizationPts_,
       fixedFields, fixedFieldValues);
}




void EquationSet::init(const Expr& eqns, 
                       const Expr& bcs,
                       const Expr& vars, 
                       const Expr& varLinearizationPts,
                       const Expr& unks,
                       const Expr& unkLinearizationPts,
                       const Expr& fixedFields,
                       const Expr& fixedFieldValues)
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

  
  /* 
   * See whether the variational functions are TestFunction objects
   * (as in a problem where we've already taken variations, or in 
   * a Galerkin-like formulation of a non-variational problem) 
   * or UnknownFunction objects, as in a variational problem. 
   */
  bool varsAreTestFunctions = false;
  for (int i=0; i<vars.size(); i++)
    {
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(vars[i].ptr().get());
      TEST_FOR_EXCEPTION(t == 0 && varsAreTestFunctions==true,
                         RuntimeError,
                         "variational function list " << vars
                         << " contains a mix of test and non-test functions");
      TEST_FOR_EXCEPTION(t != 0 && i>0 && varsAreTestFunctions==false,
                         RuntimeError,
                         "variational function list " << vars
                         << " contains a mix of test and non-test functions");
      if (t!=0) varsAreTestFunctions=true;
    }

  TEST_FOR_EXCEPTION(varsAreTestFunctions == true
                     && isVariationalProblem_,
                     RuntimeError,
                     "test functions given to a variational equation set");

  TEST_FOR_EXCEPTION(varsAreTestFunctions == false
                     && !isVariationalProblem_,
                     RuntimeError,
                     "variational functions are unknowns, but equation "
                     "set is not variational");
  
  /* map each var and unknown function's ID numbers to its
   * position in the input function lists */
  Set<int> varFuncSet;
  Set<int> unkFuncSet;
  for (int i=0; i<vars.size(); i++)
    {
      const FuncElementBase* t 
        = dynamic_cast<const FuncElementBase*>(vars[i].ptr().get());
      int fid = t->funcID();
      varFuncSet.put(fid);
      varIDToReducedIDMap_.put(fid, i);
    }

  TEST_FOR_EXCEPTION(unks.size() == 0 && !isFunctionalCalculator_,
                     InternalError,
                     "no unks passed to an equation set that is not "
                     "a gradient calculator");

  
  for (int i=0; i<unks.size(); i++)
    {
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(unks[i].ptr().get());
      TEST_FOR_EXCEPTION(u==0, RuntimeError, 
                         "EquationSet ctor input unk function "
                         << unks[i] 
                         << " does not appear to be a unk function");
      int fid = u->funcID();
      unkFuncSet.put(fid);
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
                               EvalContext::nextID(),
                               EvalContext::nextID(),
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
          varsOnRegions_.put(reg, integralSum->funcsOnRegion(d, varFuncSet));
          unksOnRegions_.put(reg, integralSum->funcsOnRegion(d, unkFuncSet));
        }
      else
        {
          Set<int>& t = varsOnRegions_.get(reg);
          t.merge(integralSum->funcsOnRegion(d, varFuncSet));
          Set<int>& u = unksOnRegions_.get(reg);
          u.merge(integralSum->funcsOnRegion(d, unkFuncSet));
        }
      for (int t=0; t<integralSum->numTerms(d); t++)
        {
          RegionQuadCombo rqc(reg.ptr(), integralSum->quad(d,t));
          Expr term = integralSum->expr(d,t);
          rqcSet.put(rqc);
          regionQuadComboExprs_.put(rqc, term);
          /* prepare calculation of both stiffness matrix and load vector */
          if (compTypes_.contains(MatrixAndVector))
            {
              EvalContext context(rqc, 2, contextID[0]);
              DerivSet nonzeros;
              
              if (isVariationalProblem_)
                {
                  nonzeros = SymbPreprocessor
                    ::setupVariations(term, 
                                      vars, varLinearizationPts,
                                      unks, unkLinearizationPts,
                                      fixedFields, fixedFieldValues,
                                      context);
                }
              else
                {
                  nonzeros = SymbPreprocessor::setupExpr(term, vars, unks, 
                                                         unkLinearizationPts,
                                                         context);
                }
              addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
                               nonzeros, false);
              rqcToContext_[MatrixAndVector].put(rqc, context);
              regionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
                                                                 nonzeros);
            }
          /* prepare calculation of load vector only */
          if (compTypes_.contains(VectorOnly))
            {
              EvalContext context(rqc, 1, contextID[1]);
              DerivSet nonzeros;
              if (isVariationalProblem_)
                {
                  nonzeros = SymbPreprocessor
                    ::setupVariations(term, vars, varLinearizationPts,
                                      unks, unkLinearizationPts,
                                      fixedFields, fixedFieldValues,
                                      context);
                }
              else
                {
                  nonzeros = SymbPreprocessor::setupExpr(term, vars, unks, 
                                                         unkLinearizationPts,
                                                         context);
                }
              rqcToContext_[VectorOnly].put(rqc, context);
              regionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
            }
          /* prepare calculation of functional value only */
          if (compTypes_.contains(FunctionalOnly))
            {
              EvalContext context(rqc, 0, contextID[2]);
              DerivSet nonzeros;
              if (compTypes_.contains(FunctionalAndGradient))
                {
                  nonzeros = SymbPreprocessor
                    ::setupGradient(term, 
                                    vars, varLinearizationPts,
                                    fixedFields, fixedFieldValues,
                                    context);
                }
              else
                {
                  nonzeros = SymbPreprocessor
                    ::setupFunctional(term, 
                                      fixedFields, fixedFieldValues,
                                      context);
                }
              rqcToContext_[FunctionalOnly].put(rqc, context);
              regionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
            }
          /* prepare calculation of functional value and gradient */
          if (compTypes_.contains(FunctionalAndGradient))
            {
              EvalContext context(rqc, 1, contextID[3]);
              DerivSet nonzeros;
              nonzeros = SymbPreprocessor
                ::setupGradient(term, 
                                vars, varLinearizationPts,
                                fixedFields, fixedFieldValues,
                                context);
              rqcToContext_[FunctionalAndGradient].put(rqc, context);
              regionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
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

          if (!regionSet.contains(reg)) 
            {
              regionSet.put(reg);
              varsOnRegions_.put(reg, bcSum->funcsOnRegion(d, varFuncSet));
              unksOnRegions_.put(reg, bcSum->funcsOnRegion(d, unkFuncSet));
              bcVarsOnRegions_.put(reg, bcSum->funcsOnRegion(d, varFuncSet));
              bcUnksOnRegions_.put(reg, bcSum->funcsOnRegion(d, unkFuncSet));
            }
          else
            {
              if (!bcVarsOnRegions_.containsKey(reg))
                {
                  bcVarsOnRegions_.put(reg, bcSum->funcsOnRegion(d, varFuncSet));
                }
              if (!bcUnksOnRegions_.containsKey(reg))
                {
                  bcUnksOnRegions_.put(reg, bcSum->funcsOnRegion(d, unkFuncSet));
                }
              Set<int>& t = varsOnRegions_.get(reg);
              t.merge(bcSum->funcsOnRegion(d, varFuncSet));
              Set<int>& u = unksOnRegions_.get(reg);
              u.merge(bcSum->funcsOnRegion(d, unkFuncSet));
            }

          for (Map<OrderedHandle<CellFilterStub>, Set<int> >::const_iterator 
                 iter=bcVarsOnRegions_.begin(); iter != bcVarsOnRegions_.end();
               iter ++)
            {
              //              cerr << "region=" << iter->first << endl;
              //              cerr << "id = " << iter->first.ptr()->id() << endl;
            }
          //          cerr << "checking map integrity " << endl;
          TEST_FOR_EXCEPTION(!bcVarsOnRegions_.containsKey(reg),
                             InternalError,
                             "Bug: region " << reg << " should appear in "
                             "BC list" << bcVarsOnRegions_);
          for (int t=0; t<bcSum->numTerms(d); t++)
            {
              RegionQuadCombo rqc(reg.ptr(), bcSum->quad(d,t));
              Expr term = bcSum->expr(d,t);
              rqcBCSet.put(rqc);
              bcRegionQuadComboExprs_.put(rqc, bcSum->expr(d,t)); 

              
              /* prepare calculation of both stiffness matrix and load vector */
              if (compTypes_.contains(MatrixAndVector))
                {
                  EvalContext context(rqc, 2, contextID[0]);
                  DerivSet nonzeros;
              
                  if (isVariationalProblem_)
                    {
                      nonzeros = SymbPreprocessor
                        ::setupVariations(term, 
                                          vars, varLinearizationPts,
                                          unks, unkLinearizationPts,
                                          fixedFields, fixedFieldValues,
                                          context);
                    }
                  else
                    {
                      nonzeros = SymbPreprocessor::setupExpr(term, vars, unks, 
                                                             unkLinearizationPts,
                                                             context);
                    }
                  addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
                                   nonzeros, true);
                  bcRqcToContext_[MatrixAndVector].put(rqc, context);
                  bcRegionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
                                                                       nonzeros);
                }
              /* prepare calculation of load vector only */
              if (compTypes_.contains(VectorOnly))
                {
                  EvalContext context(rqc, 1, contextID[1]);
                  DerivSet nonzeros;
                  if (isVariationalProblem_)
                    {
                      nonzeros = SymbPreprocessor
                        ::setupVariations(term, vars, varLinearizationPts,
                                          unks, unkLinearizationPts,
                                          fixedFields, fixedFieldValues,
                                          context);
                    }
                  else
                    {
                      nonzeros = SymbPreprocessor::setupExpr(term, vars, unks, 
                                                             unkLinearizationPts,
                                                             context);
                    }
                  bcRqcToContext_[VectorOnly].put(rqc, context);
                  bcRegionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
                }
              /* prepare calculation of functional value only */
              if (compTypes_.contains(FunctionalOnly))
                {
                  EvalContext context(rqc, 0, contextID[2]);
                  DerivSet nonzeros;
                  if (compTypes_.contains(FunctionalAndGradient))
                    {
                      nonzeros = SymbPreprocessor
                        ::setupGradient(term, 
                                        vars, varLinearizationPts,
                                        fixedFields, fixedFieldValues,
                                        context);
                    }
                  else
                    {
                      nonzeros = SymbPreprocessor
                        ::setupFunctional(term, 
                                          fixedFields, fixedFieldValues,
                                          context);
                    }
                  bcRqcToContext_[FunctionalOnly].put(rqc, context);
                  bcRegionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
                }
              /* prepare calculation of functional value and gradient */
              if (compTypes_.contains(FunctionalAndGradient))
                {
                  EvalContext context(rqc, 1, contextID[3]);
                  DerivSet nonzeros;
                  nonzeros = SymbPreprocessor
                    ::setupGradient(term, 
                                    vars, varLinearizationPts,
                                    fixedFields, fixedFieldValues,
                                    context);
                  bcRqcToContext_[FunctionalAndGradient].put(rqc, context);
                  bcRegionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
                }
            }
        }
    }

  

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Vars appearing on each region: " << endl << varsOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each region: " << endl << unksOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Vars appearing on each BC region: " 
               << endl << bcVarsOnRegions_);

  SUNDANCE_OUT(verbosity() > VerbSilent,
               "Unks appearing on each BC region: " 
               << endl << bcUnksOnRegions_);


  /* convert sets to arrays */
  regions_ = regionSet.elements();
  regionQuadCombos_ = rqcSet.elements();
  bcRegionQuadCombos_ = rqcBCSet.elements();


}


void EquationSet
::addToVarUnkPairs(const OrderedHandle<CellFilterStub>& domain,
                   const Set<int>& vars,
                   const Set<int>& unks,
                   const DerivSet& nonzeros, 
                   bool isBC)
{
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "finding var-unk pairs "
               "for domain " << domain);
  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "isBC=" << isBC);
  
  RefCountPtr<Set<OrderedPair<int, int> > > funcPairs;
  Map<OrderedHandle<CellFilterStub>, RefCountPtr<Set<OrderedPair<int, int> > > >* theMap;

  if (isBC) 
    {
      theMap = &(bcVarUnkPairsOnRegions_);
    }
  else 
    {
      theMap = &(varUnkPairsOnRegions_);
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
      
      int varID = -1;
      int unkID = -1;

      MultipleDeriv::const_iterator j;
      for (j=md.begin(); j != md.end(); j++)
        {
          const Deriv& d = *j;
          const FunctionalDeriv* fd = d.funcDeriv();
          TEST_FOR_EXCEPTION(fd==0, InternalError, "non-functional deriv "
                             << d << " detected in EquationSet::"
                             "addToVarUnkPairs()");
          if (unks.contains(fd->funcID()))
            {
              unkID = fd->funcID();
              continue;
            }
          if (vars.contains(fd->funcID()))
            {
              varID = fd->funcID();
              continue;
            }
        }
      TEST_FOR_EXCEPTION(unkID==-1 || varID==-1, InternalError,
                         "multiple derivative " << md << " does not "
                         "appear to have exactly one var and unknown each");
      funcPairs->put(OrderedPair<int, int>(varID, unkID));
    }

  SUNDANCE_OUT(verbosity() > VerbMedium, tab << "found " << *funcPairs);
  
}

const RefCountPtr<Set<OrderedPair<int, int> > >& EquationSet::
bcVarUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
{
  TEST_FOR_EXCEPTION(!bcVarUnkPairsOnRegions_.containsKey(domain),
                     InternalError,
                     "equation set does not have a var-unk pair list for "
                     "bc region " << domain);
  const RefCountPtr<Set<OrderedPair<int, int> > >& rtn 
    = bcVarUnkPairsOnRegions_.get(domain);

  TEST_FOR_EXCEPTION(rtn.get()==0, InternalError, 
                     "null var-unk pair list for BC region " << domain);
  return rtn;
}

bool EquationSet::isBCRegion(int d) const
{
  return bcVarsOnRegions_.containsKey(regions_[d]);
}


EvalContext EquationSet::rqcToContext(ComputationType compType, 
                                      const RegionQuadCombo& r) const 
{
  TEST_FOR_EXCEPTION(!rqcToContext_.containsKey(compType),
                     InternalError,
                     "EquationSet::rqcToContext() did not find key " 
                     << compType);
  return rqcToContext_.get(compType).get(r);
}

EvalContext EquationSet::bcRqcToContext(ComputationType compType, 
                                        const RegionQuadCombo& r) const 
{
  TEST_FOR_EXCEPTION(!bcRqcToContext_.containsKey(compType),
                     InternalError,
                     "EquationSet::bcRqcToContext() did not find key " 
                     << compType);
  return bcRqcToContext_.get(compType).get(r);
}

const DerivSet& EquationSet::nonzeroFunctionalDerivs(ComputationType compType,
                                                     const RegionQuadCombo& r) const
{
  TEST_FOR_EXCEPTION(!regionQuadComboNonzeroDerivs_.containsKey(compType),
                     InternalError,
                     "EquationSet:nonzeroFunctionalDerivs() did not find key " 
                     << compType);
  return regionQuadComboNonzeroDerivs_.get(compType).get(r);
}

const DerivSet& EquationSet::nonzeroBCFunctionalDerivs(ComputationType compType,
                                                       const RegionQuadCombo& r) const
{
  TEST_FOR_EXCEPTION(!bcRegionQuadComboNonzeroDerivs_.containsKey(compType),
                     InternalError,
                     "EquationSet:nonzeroBCFunctionalDerivs() did not find key " 
                     << compType);
  return bcRegionQuadComboNonzeroDerivs_.get(compType).get(r);
}
