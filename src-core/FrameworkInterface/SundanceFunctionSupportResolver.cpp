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

#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceTabs.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace std;

FunctionSupportResolver::FunctionSupportResolver(
  const Expr& eqns,
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& unks,
  const Expr& unkParams,
  const Expr& params,
  const Array<Expr>& fixedFields,
  bool isVariational,
  int verb)
  : 
    eqns_(eqns),
    bcs_(bcs),
    integralSum_(dynamic_cast<const SumOfIntegrals*>(eqns.ptr().get())),
    bcSum_(0),
    varFuncSet_(),
    unkFuncSet_(),
    unkParamSet_(),
    regions_(),
    regionToIndexMap_(),
    varsOnRegions_(),
    unksOnRegions_(),
    bcVarsOnRegions_(),
    bcUnksOnRegions_(),
    testToRegionsMap_(),
    unkToRegionsMap_(),
    varFuncData_(),
    unkFuncData_(),
    varFuncs_(flattenSpectral(vars)),
    unkFuncs_(flattenSpectral(unks)),
    fixedFields_(flattenSpectral(fixedFields)),
    unkParams_(unkParams),
    params_(params),
    varIDToReducedIDMap_(),
    unkIDToReducedIDMap_(),
    unkParamIDToReducedUnkParamIDMap_(),
    varIDToBlockMap_(),
    unkIDToBlockMap_(),
    unreducedVarID_(),
    unreducedUnkID_(),
    unreducedUnkParamID_(),
    isVariationalProblem_(isVariational)
{
  Tabs tab0;
  Tabs tab1;

  /* begin with a sanity check to ensure that the input equation set 
   * exists and is integral form */

  TEST_FOR_EXCEPTION(eqns.ptr().get()==0, RuntimeError,
    "FunctionSupportResolver ctor detected empty equation set input");

  TEST_FOR_EXCEPTION(integralSum_==0, RuntimeError,
    "FunctionSupportResolver ctor detected an input equation set that "
    "is not in integral form");

  bool hasBCs = false;
  if (bcs.ptr().get() != 0)
  {
    bcSum_ = dynamic_cast<const SumOfBCs*>(bcs.ptr().get());
    TEST_FOR_EXCEPTION(bcSum_==0, RuntimeError,
      "FunctionSupport ctor detected an input Essential "
      "BC that is not an EssentialBC object. "
      "Object is " << bcs);
    hasBCs = true;
  }

  /* upgrade base verbosity level if one of the terms is being watched */
  if (integralSum_->hasWatchedTerm() || (hasBCs && bcSum_->hasWatchedTerm()))
  {
    verb = max(verb, 1);
  }
  SUNDANCE_BANNER1(verb, tab0, "FunctionSupportResolver setup");

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


  /*
   * Now we collect arrays of the function data (e.g., basis functions).
   * Several function components may share common data. For example, the
   * components of a function discretized with Nedelec elements will
   * all point back to the same data. The symbolic engine needs to know 
   * about components, but the DOF map system only needs bases. 
   */
  varFuncData_.resize(vars.size());
  unkFuncData_.resize(unks.size());
  varIDToReducedIDMap_.resize(vars.size());
  unreducedVarID_.resize(vars.size());
  

  /* map each var and unknown function's ID numbers to its
   * position in the input function lists */
  for (unsigned int b=0; b<vars.size(); b++)
  {
    Tabs tab2;
    unreducedVarID_[b].resize(vars[b].size());
    int k=0;
    for (unsigned int i=0; i<vars[b].size(); i++, k++)
    {
      const FuncElementBase* t 
        = dynamic_cast<const FuncElementBase*>(vars[b][i].ptr().get());
      int fid = t->sharedFuncID();
      if (varFuncSet_.contains(fid)) continue;
      varFuncData_[b].append(getSharedFunctionData(t));
      varFuncSet_.put(fid);
      varIDToBlockMap_.put(fid, b);
      varIDToReducedIDMap_[b].put(fid, k);
      unreducedVarID_[b][k] = fid;
    }
    SUNDANCE_MSG2(verb, tab2 << "block=" << b << " var functions are " 
      << unreducedUnkID_[b]);
  }

  /* set up func ID maps for unks */
  unkIDToReducedIDMap_.resize(unks.size());
  unreducedUnkID_.resize(unks.size());
  for (unsigned int b=0; b<unks.size(); b++)
  {
    Tabs tab2;
    unreducedUnkID_[b].resize(unks[b].size());
    int k=0;
    for (unsigned int i=0; i<unks[b].size(); i++, k++)
    {
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(unks[b][i].ptr().get());
      TEST_FOR_EXCEPTION(u==0, RuntimeError, 
        "EquationSet ctor input unk function "
        << unks[b][i] 
        << " does not appear to be a unk function");
      int fid = u->sharedFuncID();
      if (unkFuncSet_.contains(fid)) continue;
      unkFuncData_[b].append(getSharedFunctionData(u));
      unkFuncSet_.put(fid);
      unkIDToBlockMap_.put(fid, b);
      unkIDToReducedIDMap_[b].put(fid, k);
      unreducedUnkID_[b][k] = fid;
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
    unkParamSet_.put(fid);
    unkParamIDToReducedUnkParamIDMap_.put(fid, i);
    unreducedUnkParamID_[i] = fid;
  }
  SUNDANCE_MSG2(verb, tab1 << "unk parameters are " 
    << unreducedUnkParamID_);

  Set<OrderedHandle<CellFilterStub> > regionSet;

  /* Do the non-bc eqns first */
  SUNDANCE_MSG1(verb, tab1 << "processing integral terms");
  for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
         r=integralSum_->rqcToExprMap().begin(); 
       r!=integralSum_->rqcToExprMap().end(); r++)
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


    Expr term = r->second;
    OrderedHandle<CellFilterStub> reg = rqc.domain();

    if (!regionSet.contains(reg)) 
    {
      regionSet.put(reg);
      Set<int> vf = integralSum_->funcsOnRegion(reg, varFuncSet_);
      Set<int> uf = integralSum_->funcsOnRegion(reg, unkFuncSet_);
      varsOnRegions_.put(reg, vf);
      unksOnRegions_.put(reg, uf);
    }
    else
    {
      Set<int>& t = varsOnRegions_.get(reg);
      t.merge(integralSum_->funcsOnRegion(reg, varFuncSet_));
      Set<int>& u = unksOnRegions_.get(reg);
      u.merge(integralSum_->funcsOnRegion(reg, unkFuncSet_));
    }
  }
  

  /* now do the BCs */
  if (hasBCs)
  {
    /* functions found in the BCs both in the overall lists and 
     * also in the bc-specific lists */
    for (SundanceUtils::Map<RegionQuadCombo, Expr>::const_iterator 
           r=bcSum_->rqcToExprMap().begin(); 
         r!=bcSum_->rqcToExprMap().end(); r++)
    {
      Tabs tab15;
      RegionQuadCombo rqc = r->first;
      int rqcVerb = verb;
      if (rqc.watch().isActive()) 
      {
        rqcVerb=5;
        SUNDANCE_MSG1(verb, tab15 << "processing RQC = " << rqc);
      }

      Expr term = r->second;
      OrderedHandle<CellFilterStub> reg = rqc.domain();

      if (!regionSet.contains(reg)) 
      {
        regionSet.put(reg);
        varsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        unksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
        bcVarsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        bcUnksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
      }
      else
      {
        if (!bcVarsOnRegions_.containsKey(reg))
        {
          bcVarsOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, varFuncSet_));
        }
        if (!bcUnksOnRegions_.containsKey(reg))
        {
          bcUnksOnRegions_.put(reg, bcSum_->funcsOnRegion(reg, unkFuncSet_));
        }
        Set<int>& t = varsOnRegions_.get(reg);
        t.merge(bcSum_->funcsOnRegion(reg, varFuncSet_));
        Set<int>& u = unksOnRegions_.get(reg);
        u.merge(bcSum_->funcsOnRegion(reg, unkFuncSet_));
      }
    }
  }

  /* convert sets to arrays */
  regions_ = regionSet.elements();

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




Array<Expr> FunctionSupportResolver::flattenSpectral(const Array<Expr>& expr) const
{
  Array<Expr> rtn(expr.size());
  for (unsigned int i=0; i<expr.size(); i++)
  {
    const Expr& e = expr[i];
    rtn[i] = flattenSpectral(e);
  }
  return rtn;
}

Expr FunctionSupportResolver::flattenSpectral(const Expr& expr) const
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

int FunctionSupportResolver::reducedVarID(int varID) const 
{
  TEST_FOR_EXCEPTION(!hasVarID(varID), RuntimeError, 
    "varID " << varID << " not found in equation set");

  int b = blockForVarID(varID);
  return varIDToReducedIDMap_[b].get(varID);
}

int FunctionSupportResolver::reducedUnkID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!hasUnkID(unkID), RuntimeError, 
    "unkID " << unkID << " not found in equation set");

  int b = blockForUnkID(unkID);
  return unkIDToReducedIDMap_[b].get(unkID);
}


int FunctionSupportResolver::reducedUnkParamID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!hasUnkParamID(unkID), RuntimeError, 
    "unkParamID " << unkID << " not found in equation set");

  return unkParamIDToReducedUnkParamIDMap_.get(unkID);
}

int FunctionSupportResolver::blockForVarID(int varID) const 
{
  TEST_FOR_EXCEPTION(!varIDToBlockMap_.containsKey(varID), RuntimeError,
    "key " << varID << " not found in map "
    << varIDToBlockMap_);
  return varIDToBlockMap_.get(varID);
}

int FunctionSupportResolver::blockForUnkID(int unkID) const 
{
  TEST_FOR_EXCEPTION(!unkIDToBlockMap_.containsKey(unkID), RuntimeError,
    "key " << unkID << " not found in map "
    << unkIDToBlockMap_);
  return unkIDToBlockMap_.get(unkID);
}

const Set<OrderedHandle<CellFilterStub> >&  FunctionSupportResolver::regionsForTestFunc(int testID) const
{
  TEST_FOR_EXCEPTION(!testToRegionsMap_.containsKey(testID), RuntimeError,
    "key " << testID << " not found in map "
    << testToRegionsMap_);
  return testToRegionsMap_.get(testID);
}

const Set<OrderedHandle<CellFilterStub> >&  FunctionSupportResolver::regionsForUnkFunc(int unkID) const
{
  TEST_FOR_EXCEPTION(!unkToRegionsMap_.containsKey(unkID), RuntimeError,
    "key " << unkID << " not found in map "
    << testToRegionsMap_);
  return unkToRegionsMap_.get(unkID);
}

int FunctionSupportResolver::indexForRegion(const OrderedHandle<CellFilterStub>& region) const
{
  TEST_FOR_EXCEPTION(!regionToIndexMap_.containsKey(region), RuntimeError,
    "key " << region << " not found in map "
    << regionToIndexMap_);
  return regionToIndexMap_.get(region);
}

bool FunctionSupportResolver::hasBCs() const 
{
  return bcSum_ != 0;
}


namespace SundanceCore
{
namespace Internal
{

RefCountPtr<const CommonFuncDataStub> getSharedFunctionData(const FuncElementBase* f)
{
  const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(f);
  const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(f);

  if (u != 0) 
  {
    return rcp_dynamic_cast<const CommonFuncDataStub>(u->commonData());
  }
  if (t != 0)
  {
    return rcp_dynamic_cast<const CommonFuncDataStub>(t->commonData());
  }
  
  TEST_FOR_EXCEPTION( true, InternalError, 
    "unrecognized function type: " << typeid(*f).name());
  return u->commonData(); // -Wall, will never be called;
}
}
}

