/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceEquationSet.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownParameterElement.hpp"
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
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fields,
  const Array<Expr>& fieldValues,
    const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<EquationSet>("Equation Set", verbParams),
    regions_(),
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
    testToRegionsMap_(),
    unkToRegionsMap_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    varFuncs_(),
    unkFuncs_(),
    unkLinearizationPts_(),
    unkParams_(),
    unkParamEvalPts_(),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true),
    isSensitivityProblem_(false)
{
  Array<Expr> unks;
  Array<Expr> unkEvalPt;
  Array<Expr> vars;
  Array<Expr> varEvalPt;
  
  compTypes_.put(FunctionalOnly);

  rqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(FunctionalOnly,
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalOnly,
    Map<RegionQuadCombo, DerivSet>());


  init(eqns, bcs, vars, varEvalPt,
    unks, unkEvalPt,
    unkParams_, unkParamEvalPts_,
    params, paramEvalPts,
    fields, fieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars, 
  const Array<Expr>& unks,
  const Array<Expr>& unkLinearizationPts,
  const Expr& unkParams,
  const Expr& unkParamEvalPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues,
  const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<EquationSet>("Equation Set", verbParams),
    regions_(),
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
    testToRegionsMap_(),
    unkToRegionsMap_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    varFuncs_(flattenSpectral(vars)),
    unkFuncs_(flattenSpectral(unks)),
    unkLinearizationPts_(flattenSpectral(unkLinearizationPts)),
    unkParams_(unkParams),
    unkParamEvalPts_(unkParamEvalPts),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(false),
    isFunctionalCalculator_(false),
    isSensitivityProblem_(unkParams.size() > 0)
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

  if (unkParams.size() > 0) 
  {
    compTypes_.put(Sensitivities);
    rqcToContext_.put(Sensitivities, Map<RegionQuadCombo, EvalContext>());
    bcRqcToContext_.put(Sensitivities, Map<RegionQuadCombo, EvalContext>());
    regionQuadComboNonzeroDerivs_.put(Sensitivities, 
      Map<RegionQuadCombo, DerivSet>());
    bcRegionQuadComboNonzeroDerivs_.put(Sensitivities, 
      Map<RegionQuadCombo, DerivSet>());
  }

  Array<Expr> zero;
  init(eqns, bcs, flattenSpectral(vars), zero,
    flattenSpectral(unks), flattenSpectral(unkLinearizationPts),
    unkParams_, unkParamEvalPts_,
    params, paramEvalPts,
    fixedFields, fixedFieldValues);
}


EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& varLinearizationPts, 
  const Array<Expr>& unks,
  const Array<Expr>& unkLinearizationPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues,
    const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<EquationSet>("Equation Set", verbParams), 
    regions_(),
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
    testToRegionsMap_(),
    unkToRegionsMap_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    varFuncs_(flattenSpectral(vars)),
    unkFuncs_(flattenSpectral(unks)),
    unkLinearizationPts_(flattenSpectral(unkLinearizationPts)),
    unkParams_(),
    unkParamEvalPts_(),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(false),
    isSensitivityProblem_(false)
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

  init(eqns, bcs, flattenSpectral(vars), flattenSpectral(varLinearizationPts), 
    flattenSpectral(unks), flattenSpectral(unkLinearizationPts),
    unkParams_, unkParamEvalPts_,
    params, paramEvalPts,
    fixedFields, fixedFieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& varLinearizationPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues,
    const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<EquationSet>("Equation Set", verbParams),
    regions_(),
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
    testToRegionsMap_(),
    unkToRegionsMap_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    varFuncs_(vars),
    unkFuncs_(),
    unkLinearizationPts_(),
    unkParams_(),
    unkParamEvalPts_(),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true),
    isSensitivityProblem_(false)
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

  init(eqns, bcs, flattenSpectral(vars), flattenSpectral(varLinearizationPts), 
    flattenSpectral(unkFuncs_), flattenSpectral(unkLinearizationPts_),
    unkParams_, unkParamEvalPts_,   
    params, paramEvalPts, 
    fixedFields, fixedFieldValues);
}




void EquationSet::init(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars, 
  const Array<Expr>& varLinearizationPts,
  const Array<Expr>& unks,
  const Array<Expr>& unkLinearizationPts,
  const Expr& unkParams,
  const Expr& unkParamEvalPts, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues)
{
  Tabs tab0;
  Tabs tab1;

  /* begin with a sanity check to ensure that the input equation set 
   * exists and is integral form */
  const SumOfIntegrals* integralSum
    = dynamic_cast<const SumOfIntegrals*>(eqns.ptr().get());

  TEST_FOR_EXCEPTION(eqns.ptr().get()==0, RuntimeError,
    "EquationSet ctor detected empty equation set input");

  TEST_FOR_EXCEPTION(integralSum==0, RuntimeError,
    "EquationSet ctor detected an input equation set that "
    "is not in integral form");

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
 
  int verb = verbLevel("setup");
  /* upgrade base verbosity level if one of the terms is being watched */
  if (integralSum->hasWatchedTerm() || (hasBCs && bcSum->hasWatchedTerm()))
  {
    verb = max(verb, 1);
  }
  SUNDANCE_BANNER1(verb, tab0, "EquationSet setup");

  if (hasBCs)
  {
    SUNDANCE_MSG1(verb, tab1 << "Problem has EssentialBCs");
  }
  else
  {
    SUNDANCE_MSG1(verb, tab1 << "Problem has no EssentialBCs");
  }

  /* 
   * See whether the variational functions are TestFunction objects
   * (as in a problem where we've already taken variations, or in 
   * a Galerkin-like formulation of a non-variational problem) 
   * or UnknownFunction objects, as in a variational problem. 
   */
  bool varsAreTestFunctions = false;
  for (unsigned int b=0; b<vars.size(); b++)
  {
    for (unsigned int i=0; i<vars[b].size(); i++)
    {
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(vars[b][i].ptr().get());
      TEST_FOR_EXCEPTION(t == 0 && varsAreTestFunctions==true,
        RuntimeError,
        "variational function list " << vars
        << " contains a mix of test and "
        "non-test functions");
      TEST_FOR_EXCEPTION(t != 0 && i>0 && varsAreTestFunctions==false,
        RuntimeError,
        "variational function list " << vars
        << " contains a mix of test and "
        "non-test functions");
      if (t!=0) varsAreTestFunctions=true;
    }
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

  if (isVariationalProblem_)
  {
    SUNDANCE_MSG1(verb, tab1 << "is variational problem");
  }
  else
  {
    SUNDANCE_MSG1(verb, tab1 << "is not in variational form");
  }
    
  
  /* map each var and unknown function's ID numbers to its
   * position in the input function lists */
  Set<int> varFuncSet;
  Set<int> unkFuncSet;
  Set<int> unkParamSet;
  varIDToReducedIDMap_.resize(vars.size());
  unreducedVarID_.resize(vars.size());
  for (unsigned int b=0; b<vars.size(); b++)
  {
    Tabs tab2;
    unreducedVarID_[b].resize(vars[b].size());
    for (unsigned int i=0; i<vars[b].size(); i++)
    {
      const FuncElementBase* t 
        = dynamic_cast<const FuncElementBase*>(vars[b][i].ptr().get());
      int fid = t->funcComponentID();
      varFuncSet.put(fid);
      varIDToBlockMap_.put(fid, b);
      varIDToReducedIDMap_[b].put(fid, i);
      unreducedVarID_[b][i] = fid;
    }
    SUNDANCE_MSG2(verb, tab2 << "block=" << b << " var functions are " 
      << unreducedUnkID_[b]);
  }

  TEST_FOR_EXCEPTION(unks.size() == 0 && !isFunctionalCalculator_,
    InternalError,
    "no unks passed to an equation set that is not "
    "a gradient calculator");

  /* set up func ID maps for unks */
  unkIDToReducedIDMap_.resize(unks.size());
  unreducedUnkID_.resize(unks.size());
  for (unsigned int b=0; b<unks.size(); b++)
  {
    Tabs tab2;
    unreducedUnkID_[b].resize(unks[b].size());
    for (unsigned int i=0; i<unks[b].size(); i++)
    {
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(unks[b][i].ptr().get());
      TEST_FOR_EXCEPTION(u==0, RuntimeError, 
        "EquationSet ctor input unk function "
        << unks[b][i] 
        << " does not appear to be a unk function");
      int fid = u->funcComponentID();
      unkFuncSet.put(fid);
      unkIDToBlockMap_.put(fid, b);
      unkIDToReducedIDMap_[b].put(fid, i);
      unreducedUnkID_[b][i] = fid;
    }
    SUNDANCE_MSG2(verb, tab2 << "block=" << b << " unk functions are " 
      << unreducedUnkID_[b]);
  }


  
  
  /* set up func ID maps for unk parameters */
  unreducedUnkParamID_.resize(unkParams.size());
  for (unsigned int i=0; i<unkParams.size(); i++)
  {
    const UnknownParameterElement* u 
      = dynamic_cast<const UnknownParameterElement*>(unkParams[i].ptr().get());
    TEST_FOR_EXCEPTION(u==0, RuntimeError, 
      "EquationSet ctor input unk parameter "
      << unkParams[i] 
      << " does not appear to be a unk parameter");
    int fid = u->funcComponentID();
    unkParamSet.put(fid);
    unkParamIDToReducedUnkParamIDMap_.put(fid, i);
    unreducedUnkParamID_[i] = fid;
  }
  SUNDANCE_MSG2(verb, tab1 << "unk parameters are " 
    << unreducedUnkParamID_);


  Set<OrderedHandle<CellFilterStub> > regionSet;
  Set<RegionQuadCombo> rqcSet;
  Set<RegionQuadCombo> rqcBCSet;

  Array<int> contextID = tuple(EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID());

  /* initialize list of regions to skip */
  for (Set<ComputationType>::const_iterator 
         i=compTypes_.begin(); i!=compTypes_.end(); i++)
  {
    rqcToSkip_[*i] = Set<RegionQuadCombo>();
    bcRqcToSkip_[*i] = Set<RegionQuadCombo>();
  }
  /* Now compile a list of all regions appearing in either the eqns or
   * the BCs */

  /* Do the non-bc eqns first */
  SUNDANCE_MSG1(verb, tab1 << "processing integral terms");
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         r=integralSum->rqcToExprMap().begin(); 
       r!=integralSum->rqcToExprMap().end(); r++)
  {
    Tabs tab15;
    Tabs tab2;
    RegionQuadCombo rqc = r->first;
    int rqcVerb = verb;
    if (rqc.watch().isActive()) 
    {
      rqcVerb=5;
      SUNDANCE_MSG1(verb, tab15 << "processing RQC = " << rqc);
    }


    rqcSet.put(rqc);
    Expr term = r->second;
    OrderedHandle<CellFilterStub> reg = rqc.domain();
    OrderedHandle<QuadratureFamilyStub> quad = rqc.quad();

    if (!regionSet.contains(reg)) 
    {
      regionSet.put(reg);
      Set<int> vf = integralSum->funcsOnRegion(reg, varFuncSet);
      Set<int> uf = integralSum->funcsOnRegion(reg, unkFuncSet);
      varsOnRegions_.put(reg, vf);
      unksOnRegions_.put(reg, uf);
    }
    else
    {
      Set<int>& t = varsOnRegions_.get(reg);
      t.merge(integralSum->funcsOnRegion(reg, varFuncSet));
      Set<int>& u = unksOnRegions_.get(reg);
      u.merge(integralSum->funcsOnRegion(reg, unkFuncSet));
    }

    regionQuadComboExprs_.put(rqc, term);

    /* prepare calculation of both stiffness matrix and load vector */
    if (compTypes_.contains(MatrixAndVector))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing matrix/vector calculation");
      Tabs tab3; 
      EvalContext context(rqc, 2, contextID[0]);
      DerivSet nonzeros;
      
      if (isVariationalProblem_)
      {
        nonzeros = SymbPreprocessor
          ::setupVariations(term, 
            toList(vars), 
            toList(varLinearizationPts),
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            fixedParams, 
            fixedParamEvalPts,
            context);
      }
      else
      {
        nonzeros = SymbPreprocessor
          ::setupFwdProblem(term, toList(vars), 
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context);
      }
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0U) 
      {
        rqcToSkip_[MatrixAndVector].put(rqc);
        continue;
      }
      addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
        nonzeros, false, rqcVerb);
      rqcToContext_[MatrixAndVector].put(rqc, context);
      regionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
        nonzeros);
    }



    /* prepare calculation of load vector only */
    if (compTypes_.contains(VectorOnly))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing vector-only calculation");
      Tabs tab3; 
      EvalContext context(rqc, 1, contextID[1]);
      DerivSet nonzeros;
      if (isVariationalProblem_)
      {
        nonzeros = SymbPreprocessor
          ::setupVariations(term, toList(vars), 
            toList(varLinearizationPts),
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            toList(fixedFields),
            toList(fixedFieldValues),
            fixedParams, 
            fixedParamEvalPts,
            context);
      }
      else
      {
        nonzeros = SymbPreprocessor
          ::setupFwdProblem(term, toList(vars), 
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context);
      }
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0U) 
      {
        rqcToSkip_[VectorOnly].put(rqc);
        continue;
      }
      rqcToContext_[VectorOnly].put(rqc, context);
      regionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
    }


    /* prepare calculation of sensitivities */
    if (compTypes_.contains(Sensitivities))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing sensitivity calculation");
      Tabs tab3;
      EvalContext context(rqc, 2, contextID[4]);
      DerivSet nonzeros;
      nonzeros = SymbPreprocessor
        ::setupSensitivities(term, toList(vars), 
          toList(unks), 
          toList(unkLinearizationPts),
          unkParams, 
          unkParamEvalPts,
          fixedParams, 
          fixedParamEvalPts,
          toList(fixedFields), 
          toList(fixedFieldValues),
          context);
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0U) 
      {
        rqcToSkip_[Sensitivities].put(rqc);
        continue;
      }
      rqcToContext_[Sensitivities].put(rqc, context);
      regionQuadComboNonzeroDerivs_[Sensitivities].put(rqc, nonzeros);
    }



    /* prepare calculation of functional value only */
    if (compTypes_.contains(FunctionalOnly))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing functional calculation");
      Tabs tab3;

      int maxOrder = 0;
      EvalContext context(rqc, maxOrder, contextID[2]);
      DerivSet nonzeros;
      Expr fields;
      Expr fieldValues;
      if (fixedFields.size() > 0)
      {
        if (vars.size() > 0)
        {
          fields = List(toList(fixedFields), toList(vars));
          fields = fields.flatten();
        }
        else
        {
          fields = toList(fixedFields);
        }
      }
      else
      {
        if (vars.size() > 0)
        {
          fields = toList(vars);
        }
      }
      if (fixedFieldValues.size() > 0)
      {
        if (varLinearizationPts.size() > 0)
        {
          fieldValues = List(toList(fixedFieldValues), 
            toList(varLinearizationPts));
          fieldValues = fieldValues.flatten();
        }
        else
        {
          fieldValues = toList(fixedFieldValues);
        }
      }
      else
      {
        if (varLinearizationPts.size() > 0)
        {
          fieldValues = toList(varLinearizationPts);
        }
      }

      nonzeros = SymbPreprocessor
        ::setupFunctional(term, 
          fixedParams, 
          fixedParamEvalPts,
          fields,
          fieldValues,
          context);
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);

      if (nonzeros.size()==0U) 
      {
        rqcToSkip_[FunctionalOnly].put(rqc);
        continue;
      }
      rqcToContext_[FunctionalOnly].put(rqc, context);
      regionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
    }
    /* prepare calculation of functional value and gradient */
    if (compTypes_.contains(FunctionalAndGradient))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing functional/gradient calculation");
      Tabs tab3;
      EvalContext context(rqc, 1, contextID[3]);
      DerivSet nonzeros;
      nonzeros = SymbPreprocessor
        ::setupGradient(term, 
          toList(vars), toList(varLinearizationPts),
          fixedParams, 
          fixedParamEvalPts,
          toList(fixedFields), toList(fixedFieldValues),
          context);

      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);

      if (nonzeros.size()==0U) 
      {
        rqcToSkip_[FunctionalAndGradient].put(rqc);
        continue;
      }
      rqcToContext_[FunctionalAndGradient].put(rqc, context);
      regionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
    }
  }
  
  /* now do the BCs */
  if (hasBCs)
  {
    /* functions found in the BCs both in the overall lists and 
     * also in the bc-specific lists */
    for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
           r=bcSum->rqcToExprMap().begin(); 
         r!=bcSum->rqcToExprMap().end(); r++)
    {
      Tabs tab15;
      RegionQuadCombo rqc = r->first;
      int rqcVerb = verb;
      if (rqc.watch().isActive()) 
      {
        rqcVerb=5;
        SUNDANCE_MSG1(verb, tab15 << "processing RQC = " << rqc);
      }

      rqcBCSet.put(rqc);
      Expr term = r->second;
      OrderedHandle<CellFilterStub> reg = rqc.domain();
      OrderedHandle<QuadratureFamilyStub> quad = rqc.quad();

      if (!regionSet.contains(reg)) 
      {
        regionSet.put(reg);
        varsOnRegions_.put(reg, bcSum->funcsOnRegion(reg, varFuncSet));
        unksOnRegions_.put(reg, bcSum->funcsOnRegion(reg, unkFuncSet));
        bcVarsOnRegions_.put(reg, bcSum->funcsOnRegion(reg, varFuncSet));
        bcUnksOnRegions_.put(reg, bcSum->funcsOnRegion(reg, unkFuncSet));
      }
      else
      {
        if (!bcVarsOnRegions_.containsKey(reg))
        {
          bcVarsOnRegions_.put(reg, bcSum->funcsOnRegion(reg, varFuncSet));
        }
        if (!bcUnksOnRegions_.containsKey(reg))
        {
          bcUnksOnRegions_.put(reg, bcSum->funcsOnRegion(reg, unkFuncSet));
        }
        Set<int>& t = varsOnRegions_.get(reg);
        t.merge(bcSum->funcsOnRegion(reg, varFuncSet));
        Set<int>& u = unksOnRegions_.get(reg);
        u.merge(bcSum->funcsOnRegion(reg, unkFuncSet));
      }

      TEST_FOR_EXCEPTION(!bcVarsOnRegions_.containsKey(reg),
        InternalError,
        "Bug: region " << reg << " should appear in "
        "BC list" << bcVarsOnRegions_);

      bcRegionQuadComboExprs_.put(rqc, term); 


              
      /* prepare calculation of both stiffness matrix and load vector */
      if (compTypes_.contains(MatrixAndVector))
      {
        Tabs tab3;
        EvalContext context(rqc, 2, contextID[0]);
        DerivSet nonzeros;
              
        if (isVariationalProblem_)
        {
          nonzeros = SymbPreprocessor
            ::setupVariations(term, 
              toList(vars), 
              toList(varLinearizationPts),
              toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              fixedParams, 
              fixedParamEvalPts,
              context);
        }
        else
        {
          nonzeros = SymbPreprocessor
            ::setupFwdProblem(term, toList(vars), toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              fixedParams, 
              fixedParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              context);
        }
        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0U) 
        {
          bcRqcToSkip_[MatrixAndVector].put(rqc);
          continue;
        }

        addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
          nonzeros, true, rqcVerb);
        bcRqcToContext_[MatrixAndVector].put(rqc, context);
        bcRegionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
          nonzeros);
      }





      /* prepare calculation of load vector only */
      if (compTypes_.contains(VectorOnly))
      {
        Tabs tab3;
        EvalContext context(rqc, 1, contextID[1]);
        DerivSet nonzeros;
        if (isVariationalProblem_)
        {
          nonzeros = SymbPreprocessor
            ::setupVariations(term, toList(vars), 
              toList(varLinearizationPts),
              toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              fixedParams, 
              fixedParamEvalPts,
              context);
        }
        else
        {
          nonzeros = SymbPreprocessor
            ::setupFwdProblem(term, toList(vars), toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              fixedParams, 
              fixedParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              context);
        }
        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0U) 
        {
          bcRqcToSkip_[VectorOnly].put(rqc);
          continue;
        }
        bcRqcToContext_[VectorOnly].put(rqc, context);
        bcRegionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
      }








      /* prepare calculation of sensitivities */
      if (compTypes_.contains(Sensitivities))
      {
        Tabs tab3;
        EvalContext context(rqc, 2, contextID[4]);
        DerivSet nonzeros;
        nonzeros = SymbPreprocessor
          ::setupSensitivities(term, toList(vars), toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0U) 
        {
          bcRqcToSkip_[Sensitivities].put(rqc);
          continue;
        }
        bcRqcToContext_[Sensitivities].put(rqc, context);
        bcRegionQuadComboNonzeroDerivs_[Sensitivities].put(rqc, nonzeros);
      }







      /* prepare calculation of functional value only */
      if (compTypes_.contains(FunctionalOnly))
      {
        Tabs tab3;
        EvalContext context(rqc, 0, contextID[2]);
        DerivSet nonzeros;
        Expr fields;
        Expr fieldValues;
        if (fixedFields.size() > 0)
        {
          if (vars.size() > 0)
          {
            fields = List(toList(fixedFields), toList(vars));
            fields = fields.flatten();
          }
          else
          {
            fields = toList(fixedFields);
          }
        }
        else
        {
          if (vars.size() > 0)
          {
            fields = toList(vars);
          }
        }
        if (fixedFieldValues.size() > 0)
        {
          if (varLinearizationPts.size() > 0)
          {
            fieldValues = List(toList(fixedFieldValues), 
              toList(varLinearizationPts));
            fieldValues = fieldValues.flatten();
          }
          else
          {
            fieldValues = toList(fixedFieldValues);
          }
        }
        else
        {
          if (varLinearizationPts.size() > 0)
          {
            fieldValues = toList(varLinearizationPts);
          }
        }
        nonzeros = SymbPreprocessor
          ::setupFunctional(term, 
            fixedParams, 
            fixedParamEvalPts,
            fields, fieldValues,
            context);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);                  
        if (nonzeros.size()==0U) 
        {
          bcRqcToSkip_[FunctionalOnly].put(rqc);
          continue;
        }
        bcRqcToContext_[FunctionalOnly].put(rqc, context);
        bcRegionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
      }




      /* prepare calculation of functional value and gradient */
      if (compTypes_.contains(FunctionalAndGradient))
      {
        Tabs tab3;
        EvalContext context(rqc, 1, contextID[3]);
        DerivSet nonzeros;
        nonzeros = SymbPreprocessor
          ::setupGradient(term, 
            toList(vars), 
            toList(varLinearizationPts),
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0U) 
        {
          bcRqcToSkip_[FunctionalAndGradient].put(rqc);
          continue;
        }
        bcRqcToContext_[FunctionalAndGradient].put(rqc, context);
        bcRegionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
      }
    }
  }

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
    "Vars appearing on each region: " << endl << varsOnRegions_);

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
    "Unks appearing on each region: " << endl << unksOnRegions_);

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
    "Vars appearing on each BC region: " 
    << endl << bcVarsOnRegions_);

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
    "Unks appearing on each BC region: " 
    << endl << bcUnksOnRegions_);


  /* convert sets to arrays */
  regions_ = regionSet.elements();
  regionQuadCombos_ = rqcSet.elements();
  bcRegionQuadCombos_ = rqcBCSet.elements();

  reducedVarsOnRegions_.resize(regions_.size());
  reducedUnksOnRegions_.resize(regions_.size());
  for (unsigned int r=0; r<regions_.size(); r++)
  {
    regionToIndexMap_.put(regions_[r], r);
    reducedVarsOnRegions_[r].resize(numVarBlocks());
    reducedUnksOnRegions_[r].resize(numUnkBlocks());
    OrderedHandle<CellFilterStub> cf = regions_[r];
    const Set<int>& v = this->varsOnRegion(r);
    const Set<int>& u = this->unksOnRegion(r);
    const Set<int>& bv = this->varsOnRegion(r);
    const Set<int>& bu = this->unksOnRegion(r);
    Set<int> vf = v;
    Set<int> uf = u;
    vf.merge(bv);
    uf.merge(bu);
    for (Set<int>::const_iterator i=vf.begin(); i!=vf.end(); i++)
    {
      int fid = *i;
      if (testToRegionsMap_.containsKey(fid))
      {
        testToRegionsMap_[fid].put(cf);
      }
      else
      {
        Set<OrderedHandle<CellFilterStub> > s;
        s.put(cf);
        testToRegionsMap_.put(fid, s);
      }
      int rv = reducedVarID(fid);
      int br = blockForVarID(fid);
      reducedVarsOnRegions_[r][br].put(rv);
    }
    for (Set<int>::const_iterator i=uf.begin(); i!=uf.end(); i++)
    {
      int fid = *i;
      if (unkToRegionsMap_.containsKey(fid))
      {
        unkToRegionsMap_[fid].put(cf);
      }
      else
      {
        Set<OrderedHandle<CellFilterStub> > s;
        s.put(cf);
        unkToRegionsMap_.put(fid, s);
      }
      int ru = reducedUnkID(fid);
      int bc = blockForUnkID(fid);
      reducedUnksOnRegions_[r][bc].put(ru);
    }
  }


}


void EquationSet
::addToVarUnkPairs(const OrderedHandle<CellFilterStub>& domain,
  const Set<int>& vars,
  const Set<int>& unks,
  const DerivSet& nonzeros, 
  bool isBC,
  int verb)
{
  Tabs tab;
  SUNDANCE_MSG2(verb, tab << "finding var-unk pairs "
    "for domain " << domain);
  SUNDANCE_MSG2(verb, tab << "isBC=" << isBC);
  
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
      
    Array<const FunctionalDeriv*> f;
    for (MultipleDeriv::const_iterator j=md.begin(); j != md.end(); j++)
    {
      const Deriv& d = *j;
      const FunctionalDeriv* fd = d.funcDeriv();
      TEST_FOR_EXCEPTION(fd==0, InternalError, "non-functional deriv "
        << d << " detected in EquationSet::"
        "addToVarUnkPairs()");
      f.append(fd);
    }
    
    bool gotIt=false;
    if (unks.contains(f[0]->funcComponentID())
      && vars.contains(f[1]->funcComponentID()))
    {
      int unkID = f[0]->funcComponentID();
      int varID = f[1]->funcComponentID();
      funcPairs->put(OrderedPair<int, int>(varID, unkID));
      gotIt=true;
    }
    if (unks.contains(f[1]->funcComponentID())
      && vars.contains(f[0]->funcComponentID()))
    {
      int unkID = f[1]->funcComponentID();
      int varID = f[0]->funcComponentID();
      funcPairs->put(OrderedPair<int, int>(varID, unkID));
      gotIt=true;
    }
    TEST_FOR_EXCEPTION(!gotIt, InternalError,
      "no valid (var,unk) pair could be extracted from "
      "derivative " << md);
  }

  SUNDANCE_MSG2(verb, tab << "found " << *funcPairs);
  
}

bool EquationSet::hasActiveWatchFlag() const 
{
  for (unsigned int i=0; i<regionQuadCombos().size(); i++)
  {
    if (regionQuadCombos()[i].watch().isActive()) return true;
  }
  for (unsigned int i=0; i<bcRegionQuadCombos().size(); i++)
  {
    if (bcRegionQuadCombos()[i].watch().isActive()) return true;
  }
  return false;
}

Array<Expr> EquationSet::flattenSpectral(const Array<Expr>& expr) const
{
  Array<Expr> rtn(expr.size());
  for (unsigned int i=0; i<expr.size(); i++)
  {
    const Expr& e = expr[i];
    rtn[i] = flattenSpectral(e);
  }
  return rtn;
}

Expr EquationSet::flattenSpectral(const Expr& expr) const
{
  Array<Expr> rtn(expr.size());
  for (unsigned int i=0; i<expr.size(); i++)
  {
    if (expr[i].size() == 1)
    {
      const SpectralExpr* se 
        = dynamic_cast<const SpectralExpr*>(expr[i][0].ptr().get());
      if (se != 0)
      {
        int nt = se->getSpectralBasis().nterms();
        Array<Expr> e(nt);
        for (int j=0; j<nt; j++)
        {
          e[j] = se->getCoeff(j);
        }
        rtn[i] = new ListExpr(e);
      }
      else
      {
        rtn[i] = expr[i];
      }
    }
    else
    {
      rtn[i] = flattenSpectral(expr[i]);
    }
  }
  Expr r = new ListExpr(rtn);
  return r.flatten();
                  
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
  TEST_FOR_EXCEPTION(!rqcToContext_.get(compType).containsKey(r),
    InternalError, 
    "EquationSet::rqcToContext(" << compType 
    << ") did not find expected key " 
    << r);
    
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


bool EquationSet::skipRqc(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEST_FOR_EXCEPTION(!rqcToSkip_.containsKey(compType),
    InternalError,
    "EquationSet::skipRqc() did not find expected key " 
    << compType);
    
  return rqcToSkip_.get(compType).contains(r);
}

bool EquationSet::skipBCRqc(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEST_FOR_EXCEPTION(!bcRqcToSkip_.containsKey(compType),
    InternalError,
    "EquationSet::skipBCRqc() did not find expected key " 
    << compType);
    
  return bcRqcToSkip_.get(compType).contains(r);
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



int EquationSet::reducedVarID(int varID) const 
{
  TEST_FOR_EXCEPTION(!hasVarID(varID), RuntimeError, 
    "varID " << varID << " not found in equation set");

  int b = blockForVarID(varID);
  return varIDToReducedIDMap_[b].get(varID);
}

int EquationSet::reducedUnkID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!hasUnkID(unkID), RuntimeError, 
    "unkID " << unkID << " not found in equation set");

  int b = blockForUnkID(unkID);
  return unkIDToReducedIDMap_[b].get(unkID);
}


int EquationSet::reducedUnkParamID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!hasUnkParamID(unkID), RuntimeError, 
    "unkParamID " << unkID << " not found in equation set");

  return unkParamIDToReducedUnkParamIDMap_.get(unkID);
}


Expr EquationSet::toList(const Array<Expr>& e)
{
  return new ListExpr(e);
}


int EquationSet::blockForVarID(int varID) const 
{
  TEST_FOR_EXCEPTION(!varIDToBlockMap_.containsKey(varID), RuntimeError,
    "key " << varID << " not found in map "
    << varIDToBlockMap_);
  return varIDToBlockMap_.get(varID);
}

int EquationSet::blockForUnkID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!unkIDToBlockMap_.containsKey(unkID), RuntimeError,
    "key " << unkID << " not found in map "
    << unkIDToBlockMap_);
  return unkIDToBlockMap_.get(unkID);
}

const Set<OrderedHandle<CellFilterStub> >&  EquationSet::regionsForTestFunc(int testID) const
{
  TEST_FOR_EXCEPTION(!testToRegionsMap_.containsKey(testID), RuntimeError,
    "key " << testID << " not found in map "
    << testToRegionsMap_);
  return testToRegionsMap_.get(testID);
}

const Set<OrderedHandle<CellFilterStub> >&  EquationSet::regionsForUnkFunc(int unkID) const
{
  TEST_FOR_EXCEPTION(!unkToRegionsMap_.containsKey(unkID), RuntimeError,
    "key " << unkID << " not found in map "
    << testToRegionsMap_);
  return unkToRegionsMap_.get(unkID);
}

int EquationSet::indexForRegion(const OrderedHandle<CellFilterStub>& region) const
{
  TEST_FOR_EXCEPTION(!regionToIndexMap_.containsKey(region), RuntimeError,
    "key " << region << " not found in map "
    << regionToIndexMap_);
  return regionToIndexMap_.get(region);
}
