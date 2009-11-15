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

#include "SundanceAssembler.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceTrivialGrouper.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_HashTable.h"
#include "SundanceIntHashSet.hpp"
#include "TSFProductVectorSpaceDecl.hpp"
#include "TSFLoadableBlockVector.hpp"
#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFBlockOperatorBaseDecl.hpp"
#include "TSFSimpleBlockOpDecl.hpp"
#include "SundanceAssemblyKernelBase.hpp"
#include "SundanceVectorAssemblyKernel.hpp"
#include "SundanceMatrixVectorAssemblyKernel.hpp"
#include "SundanceFunctionalAssemblyKernel.hpp"
#include "SundanceFunctionalGradientAssemblyKernel.hpp"
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFLinearOperatorImpl.hpp"
#include "TSFSimpleBlockOpImpl.hpp"
#endif



using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

static Time& assemblyTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("assembly"); 
  return *rtn;
}

static Time& assemblerCtorTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("assembler ctor"); 
  return *rtn;
}


static Time& configTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("matrix config"); 
  return *rtn;
}

static Time& graphBuildTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("matrix graph determination"); 
  return *rtn;
}

static Time& colSearchTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("graph column processing"); 
  return *rtn;
}



static Time& matAllocTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("matrix allocation"); 
  return *rtn;
}

static Time& matFinalizeTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("matrix graph packing"); 
  return *rtn;
}

static Time& graphFlatteningTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("tmp graph flattening"); 
  return *rtn;
}


RefCountPtr<ParameterList> Assembler::defaultVerbParams()
{
  static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Assembler"));
  static int first = true;
  if (first)
  {
    rtn->set<int>("setup", 0);
    rtn->set<int>("matrix config", 0);
    rtn->set<int>("vector config", 0);
    rtn->set("DOF Map", *DOFMapBase::defaultVerbParams());
    first = false;
  }
  return rtn;
}

Assembler
::Assembler(const Mesh& mesh, 
  const RefCountPtr<EquationSet>& eqn,
  const Array<VectorType<double> >& rowVectorType,
  const Array<VectorType<double> >& colVectorType,
  bool partitionBCs,
  const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<Assembler>("Assembler", verbParams),
    partitionBCs_(partitionBCs),
    matNeedsConfiguration_(true),
    matNeedsFinalization_(true),
    numConfiguredColumns_(0),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    externalRowSpace_(eqn->numVarBlocks()),
    externalColSpace_(eqn->numUnkBlocks()),
    privateRowSpace_(eqn->numVarBlocks()),
    privateColSpace_(eqn->numUnkBlocks()),
    bcRows_(eqn->numVarBlocks()),
    bcCols_(eqn->numUnkBlocks()),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    isInternalBdry_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(eqn->numVarBlocks()),
    isBCCol_(eqn->numUnkBlocks()),
    remoteBCCols_(eqn->numUnkBlocks()),
    lowestRow_(eqn->numVarBlocks()),
    lowestCol_(eqn->numUnkBlocks()),
    rowVecType_(rowVectorType),
    colVecType_(colVectorType),
    testIDToBlockMap_(),
    unkIDToBlockMap_(),
    converter_(eqn->numUnkBlocks())
{
  TimeMonitor timer(assemblerCtorTimer());
  init(mesh, eqn);
}

Assembler
::Assembler(const Mesh& mesh, 
  const RefCountPtr<EquationSet>& eqn,
  const ParameterList& verbParams)
  : ParameterControlledObjectWithVerbosity<Assembler>("Assembler", verbParams),
    partitionBCs_(false),
    matNeedsConfiguration_(true),
    matNeedsFinalization_(true),
    numConfiguredColumns_(0),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    externalRowSpace_(eqn->numVarBlocks()),
    externalColSpace_(eqn->numUnkBlocks()),
    privateRowSpace_(eqn->numVarBlocks()),
    privateColSpace_(eqn->numUnkBlocks()),
    bcRows_(eqn->numVarBlocks()),
    bcCols_(eqn->numUnkBlocks()),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    isInternalBdry_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(eqn->numVarBlocks()),
    isBCCol_(eqn->numUnkBlocks()),
    remoteBCCols_(eqn->numUnkBlocks()),
    lowestRow_(eqn->numVarBlocks()),
    lowestCol_(eqn->numUnkBlocks()),
    rowVecType_(),
    colVecType_(),
    testIDToBlockMap_(),
    unkIDToBlockMap_(),
    fixedParamIDToVectorNumber_(),
    converter_(eqn->numUnkBlocks())
{
  TimeMonitor timer(assemblerCtorTimer());
  init(mesh, eqn);
}

void Assembler::init(const Mesh& mesh, 
  const RefCountPtr<EquationSet>& eqn)
{
  Tabs tab0(0);

  /* Decide a verbosity level for the overall setup */
  int verb = this->verbLevel("setup");
  if (eqn->hasActiveWatchFlag()) 
  {
    verb = max(verb, 1);
  }

  SUNDANCE_BANNER1(verb, tab0, "Assembler setup");


  /* Create an integral grouper */
  RefCountPtr<GrouperBase> grouper 
    = rcp(new TrivialGrouper());


  /* Find out which types of computations this assembler will 
   * be required to carry out */
  const Set<ComputationType>& compTypes = eqn->computationTypes();


  /* Create the DOF map for the row space */
  DOFMapBuilder mapBuilder;

  if (compTypes.contains(VectorOnly) 
    || compTypes.contains(Sensitivities) 
    || compTypes.contains(FunctionalAndGradient))
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "building row spaces");
    mapBuilder = DOFMapBuilder(mesh, eqn->fsr(), partitionBCs_, 
      verbSublist("DOF Map"));

    rowMap_ = mapBuilder.rowMap();
    isBCRow_ = mapBuilder.isBCRow();
    isBCCol_ = mapBuilder.isBCCol();
    lowestRow_.resize(eqn_->numVarBlocks());
    /* create discrete space for each block */
    for (unsigned int b=0; b<eqn_->numVarBlocks(); b++) 
    {
      Tabs tab2;
      lowestRow_[b] = rowMap_[b]->lowestLocalDOF();
      SUNDANCE_MSG2(verb, tab2 << "block " << b << ": lowest row="
        << lowestRow_[b]);
      externalRowSpace_[b] = rcp(
        new DiscreteSpace(mesh, testBasisArray(mapBuilder.fsr())[b], 
          rowMap_[b], rowVecType_[b]));
      if (partitionBCs_)
      {
        privateRowSpace_[b] = rcp(
          new DiscreteSpace(mesh, testBasisArray(mapBuilder.fsr())[b], 
            rowMap_[b], isBCRow_[b], rowVecType_[b]));
      }
      else
      {
        privateRowSpace_[b] = externalRowSpace_[b];
      }
      SUNDANCE_MSG2(verb, tab2 << "block " << b << ": done forming row space");
    }
  }

  if (!eqn->isFunctionalCalculator())
  {
    Tabs tab1;
    /* Create the DOF map for the column space */
    SUNDANCE_MSG2(verb, tab1 << "building column spaces");
    colMap_ = mapBuilder.colMap();
    /* create discrete space for each block */
    for (unsigned int b=0; b<eqn_->numUnkBlocks(); b++) 
    {
      Tabs tab2;
      externalColSpace_[b] 
        = rcp(new DiscreteSpace(mesh, unkBasisArray(mapBuilder.fsr())[b], 
            colMap_[b], colVecType_[b]));
      if (partitionBCs_)
      {
        privateColSpace_[b] 
          = rcp(new DiscreteSpace(mesh, unkBasisArray(mapBuilder.fsr())[b], 
              colMap_[b], isBCCol_[b], colVecType_[b]));
        converter_[b] 
          = rcp(new PartitionedToMonolithicConverter(
                  privateColSpace_[b]->vecSpace(), 
                  isBCCol_[b], externalColSpace_[b]->vecSpace()));
      }
      else
      {
        privateColSpace_[b] = externalColSpace_[b];
      }
      SUNDANCE_MSG2(verb, tab2 << "block " << b << ": done forming col space");
    }

    /* initialize empty tables of information for each RQC in a 
     * matrix-vector calculation */
    groups_.put(MatrixAndVector, Array<Array<RCP<IntegralGroup> > >());
    rqcRequiresMaximalCofacets_.put(MatrixAndVector, 
      Array<IntegrationCellSpecifier>());
    contexts_.put(MatrixAndVector, Array<EvalContext>());
    evalExprs_.put(MatrixAndVector, Array<const EvaluatableExpr*>());

    /* create tables for vector calculation */
    groups_.put(VectorOnly, Array<Array<RCP<IntegralGroup> > >());
    rqcRequiresMaximalCofacets_.put(VectorOnly, 
      Array<IntegrationCellSpecifier>());
    contexts_.put(VectorOnly, Array<EvalContext>());
    evalExprs_.put(VectorOnly, Array<const EvaluatableExpr*>());

    if (eqn->isSensitivityCalculator())
    {
      fixedParamIDToVectorNumber_ 
        = eqn->fsr()->fixedParamIDToReducedFixedParamIDMap();

      /* create tables for sensitivity calculation */
      groups_.put(Sensitivities, Array<Array<RCP<IntegralGroup> > >());
      rqcRequiresMaximalCofacets_.put(Sensitivities, 
        Array<IntegrationCellSpecifier>());
      contexts_.put(Sensitivities, Array<EvalContext>());
      evalExprs_.put(Sensitivities, Array<const EvaluatableExpr*>());
      
    }
  }
  else
  {
    /* create tables for functional and gradient calculation */
    groups_.put(FunctionalAndGradient, Array<Array<RCP<IntegralGroup> > >());
    rqcRequiresMaximalCofacets_.put(FunctionalAndGradient, 
      Array<IntegrationCellSpecifier>());
    contexts_.put(FunctionalAndGradient, Array<EvalContext>());
    evalExprs_.put(FunctionalAndGradient, Array<const EvaluatableExpr*>());
    /* create tables for functional calculation */
    groups_.put(FunctionalOnly, Array<Array<RCP<IntegralGroup> > >());
    rqcRequiresMaximalCofacets_.put(FunctionalOnly, 
      Array<IntegrationCellSpecifier>());
    contexts_.put(FunctionalOnly, Array<EvalContext>());
    evalExprs_.put(FunctionalOnly, Array<const EvaluatableExpr*>());
  }





  /* --- We now loop over non-BC RQCs, doing initialization tasks for each */
  SUNDANCE_MSG1(verb, tab0 << endl 
    << tab0 << "=== setting up non-BC region-quadrature combinations");

  for (unsigned int r=0; r<eqn->regionQuadCombos().size(); r++)
  {
    Tabs tab1;
    Tabs tab12;
    const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];

    /* Determine the verbosity setting for this RQC */
    bool watchMe = rqc.watch().isActive();

    int rqcVerb = verb;
    int integralCtorVerb = 0;
    int integrationVerb = 0;
    int integralTransformVerb = 0;
    if (watchMe) 
    {
      rqcVerb = max(4,rqcVerb);
      integralCtorVerb = rqc.watch().param("integration setup");
      integrationVerb = rqc.watch().param("integration");
      integralTransformVerb = rqc.watch().param("integral transformation");
    }


    /* Note that I'm not an essential BC */
    rqc_.append(rqc);
    isBCRqc_.append(false);

    /* Find the expression for this RQC */
    const Expr& expr = eqn->expr(rqc);

    SUNDANCE_MSG2(rqcVerb, tab1 << endl << tab1 << "------------------------------------------------");
    SUNDANCE_MSG2(rqcVerb, tab1 << "initializing assembly for"
      << endl << tab12 << "rqc=" 
      << rqc << endl << tab12 << endl << tab12 << "------------------------------"
      << endl << tab12 << "expr = " << expr
      << endl << tab12 << "------------------------------"
      );

    
    /* Find the cell type needed for this RQC */
    int cellDim = CellFilter(rqc.domain()).dimension(mesh);
    CellType cellType = mesh.cellType(cellDim);
    CellType maxCellType = mesh.cellType(mesh.spatialDim());
    QuadratureFamily quad(rqc.quad());

    /* Detect internal boundaries. These need special handling */
    bool isInternalBdry = detectInternalBdry(cellDim, rqc.domain());
    isInternalBdry_.append(isInternalBdry);

    SUNDANCE_MSG2(rqcVerb, tab12 << "isInternalBdry=" << isInternalBdry);

    /* Do setup for each required computation type */
    bool rqcUsed = false;

    for (Set<ComputationType>::const_iterator 
           i=eqn->computationTypes().begin(); 
         i!=eqn->computationTypes().end();
         i++)
    {
      Tabs tab2;
      const ComputationType& compType = *i;
      SUNDANCE_MSG2(rqcVerb, tab12 << endl << tab12
        << "** computation type " << compType);
      /* Some RQCs may be unused in a given calculation. For example, an RQC
       * may be needed for vector calculation but not matrix-vector 
       * calculation. See if this RQC is needed for the present 
       * computation type. If not, there's nothing more to do here. */
      if (eqn->skipRqc(compType, rqc))
      {
        SUNDANCE_MSG2(rqcVerb, tab2 << "RQC not needed for computation type  " 
          << compType);
        continue;
      }

      /* If we're to this point, we know the RQC is needed for this 
       * computation type */
      rqcUsed = true;

      /* Look up a "context" object that we'll use as a key for different
       * evaluations of this expression */
      EvalContext context = eqn->rqcToContext(compType, rqc);

      SUNDANCE_MSG2(rqcVerb, tab2 << "context " << context.brief());

      /* Register the context */
      contexts_[compType].append(context);

      /* Register the expression */
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_[compType].append(ee);

      /* Get the "sparsity superset" which is a description of all 
       * derivatives that must be computed by this expression in the
       * present context */
      const RefCountPtr<SparsitySuperset>& sparsity 
        = ee->sparsitySuperset(context);
      SUNDANCE_MSG3(rqcVerb, tab2 << "sparsity pattern " << *sparsity);

      /* We're now ready to create integration groups for doing the 
       * integrals needed in this computation for the present RQC. */
      Array<RCP<IntegralGroup> > groups;
      grouper->setVerbosity(integralCtorVerb, integrationVerb, integralTransformVerb);
      grouper->findGroups(*eqn, maxCellType, mesh.spatialDim(),
        cellType, cellDim, quad, sparsity, isInternalBdry, groups);
      grouper->setVerbosity(0,0,0);
      groups_[compType].append(groups);

      /* Record whether or not integrations need to be done by reference
       * to maximal cofacets. */ 
      IntegrationCellSpecifier cellSpec 
        = whetherToUseCofacets(groups, ee, 
          cellDim==mesh_.spatialDim(), rqcVerb);
      SUNDANCE_MSG2(rqcVerb, tab2 << "integration: " << cellSpec);
      rqcRequiresMaximalCofacets_[compType].append(cellSpec);
      SUNDANCE_MSG2(rqcVerb, tab12
        << "done with computation type " << compType);
    }
    
    /* If this RQC has never been used, we've made a mistake */
    TEST_FOR_EXCEPTION(!rqcUsed, InternalError, "rqc=" << rqc 
      << " never used for any computation???");
    SUNDANCE_MSG2(rqcVerb, tab12 << "creating evaluation mediator for rqc=" 
      << rqc << endl << tab12 << "expr = " << expr);

    /* Finally, create an evaluation mediator for this RQC. The evaluation
     * mediator is the object through which symbolic objects refer to
     * mesh-dependent quantities (e.g., discrete functions) during
     * evaluation.  */
    mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
          quad)));
    SUNDANCE_MSG2(rqcVerb, tab1 
      << "done with RQC");
  }



  /* --- We now loop over BC RQCs, doing initialization tasks for each */
  SUNDANCE_MSG1(verb, tab0 << endl 
    << tab0 << "=== setting up BC region-quadrature combinations");
  
  for (unsigned int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
  {
    Tabs tab1;
    const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];

    /* Determine the verbosity setting for this RQC */
    bool watchMe = rqc.watch().isActive();
    int rqcVerb = verb;
    int integralCtorVerb = 0;
    int integrationVerb = 0;
    int integralTransformVerb = 0;
    if (watchMe) 
    {
      rqcVerb = max(4,rqcVerb);
      integralCtorVerb = rqc.watch().param("integration setup");
      integrationVerb = rqc.watch().param("integration");
      integralTransformVerb = rqc.watch().param("integral transformation");
    }

    /* Note that I am an essential BC */
    rqc_.append(rqc);
    isBCRqc_.append(true);


    /* Find the expression for this RQC */    
    const Expr& expr = eqn->bcExpr(rqc);

    SUNDANCE_MSG2(rqcVerb, tab1 << endl << tab1 
      << "------------------------------------------------");
    SUNDANCE_MSG1(rqcVerb, tab1 << "initializing assembly for BC rqc=" 
      << rqc << endl << tab1 
      << "expr = " << expr);
      
    /* Find the cell type needed for this RQC */    
    int cellDim = CellFilter(rqc.domain()).dimension(mesh);
    CellType cellType = mesh.cellType(cellDim);
    CellType maxCellType = mesh.cellType(mesh.spatialDim());
    QuadratureFamily quad(rqc.quad());

    /* Detect internal boundaries. These need special handling */
    bool isInternalBdry = detectInternalBdry(cellDim, rqc.domain());
    isInternalBdry_.append(isInternalBdry);

    /* Do setup for each required computation type */
    bool rqcUsed = false;

    for (Set<ComputationType>::const_iterator 
           i=eqn->computationTypes().begin(); 
         i!=eqn->computationTypes().end();
         i++)
    {
      Tabs tab2;
      const ComputationType& compType = *i;
      SUNDANCE_MSG2(rqcVerb, tab1 << endl << tab1 
        << "** computation type " << compType);
      
      /* Some RQCs may be unused in a given calculation. For example, an RQC
       * may be needed for vector calculation but not matrix-vector 
       * calculation. See if this RQC is needed for the present 
       * computation type. If not, there's nothing more to do here. */
      if (eqn->skipBCRqc(compType, rqc))
      {
        SUNDANCE_MSG2(rqcVerb, 
          tab2 << "this rqc not needed for computation type " << compType);
        continue;
      }

      /* If we're to this point, we know the RQC is needed for this 
       * computation type */   
      rqcUsed = true;

      /* Look up a "context" object that we'll use as a key for different
       * evaluations of this expression */
      EvalContext context = eqn->bcRqcToContext(compType, rqc);
      SUNDANCE_MSG2(rqcVerb, tab2 << "context " << context);

      
      contexts_[compType].append(context);
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_[compType].append(ee);
      const RefCountPtr<SparsitySuperset>& sparsity 
        = ee->sparsitySuperset(context);
      SUNDANCE_MSG3(rqcVerb, tab2 << "sparsity pattern " << *sparsity);

      Array<RCP<IntegralGroup> > groups;
      grouper->setVerbosity(integralCtorVerb, integrationVerb, integralTransformVerb);
      grouper->findGroups(*eqn, maxCellType, mesh.spatialDim(),
        cellType, cellDim, quad, sparsity, isInternalBdry, groups);
      grouper->setVerbosity(0,0,0);
      groups_[compType].append(groups);
      IntegrationCellSpecifier cellSpec 
        = whetherToUseCofacets(groups, ee, 
          cellDim==mesh_.spatialDim(), rqcVerb);
      SUNDANCE_MSG2(rqcVerb, tab2 << "integration: " << cellSpec);
      rqcRequiresMaximalCofacets_[compType].append(cellSpec);
      SUNDANCE_MSG2(rqcVerb, tab1
        << "done with computation type " << compType);
    }
    TEST_FOR_EXCEPTION(!rqcUsed, InternalError, "BC rqc=" << rqc 
      << " never used for any computation???");


    SUNDANCE_MSG2(rqcVerb, tab1 << "creating evaluation mediator for BC rqc=" 
      << rqc << endl << tab1 << "expr = " << expr);
    mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
          quad)));
    SUNDANCE_MSG2(rqcVerb, tab1 
      << "done with BC RQC");
  }
}

bool Assembler::detectInternalBdry(int cellDim,
  const CellFilter& filter) const
{
  int d = mesh_.spatialDim();
  if (cellDim == d-1)
  {
    CellSet cells = filter.getCells(mesh_);
    for (CellIterator c=cells.begin(); c!=cells.end(); c++)
    {
      if (mesh_.numMaxCofacets(cellDim, *c) > 1) return true;
    }      
  }
  return false;
}

IntegrationCellSpecifier Assembler::whetherToUseCofacets(
  const Array<RCP<IntegralGroup> >& groups,
  const EvaluatableExpr* ee,
  bool isMaximalCell,
  int verb) const
{
  Tabs tab;
  SUNDANCE_MSG2(verb, 
    tab << "deciding whether to use cofacet cells for some integrations");

  if (isMaximalCell)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, 
      tab1 << "cofacets not needed because cells are maximal");
    return NoTermsNeedCofacets;
  }
  
  IntegrationCellSpecifier cellSpec = SomeTermsNeedCofacets;

  bool allTermsNeedCofacets = true;
  bool noTermsNeedCofacets = true;
  for (unsigned int g=0; g<groups.size(); g++)
  {
    Tabs tab1;
    switch(groups[g]->usesMaximalCofacets()) 
    {
      case NoTermsNeedCofacets:
        allTermsNeedCofacets = false;
        SUNDANCE_MSG2(verb, 
          tab1 << "integral group " << g << " does not need cofacets");
        break;
      case AllTermsNeedCofacets:
      case SomeTermsNeedCofacets:
        noTermsNeedCofacets = false;
        SUNDANCE_MSG2(verb, 
          tab1 << "integral group " << g << " needs cofacets");
        break;
      default:
        TEST_FOR_EXCEPT(1);
    }
  } 

  if (allTermsNeedCofacets)
  {
    cellSpec = AllTermsNeedCofacets;
  }
  else if (noTermsNeedCofacets)
  {
    cellSpec = NoTermsNeedCofacets;
  }

  if (!isMaximalCell && ee->maxDiffOrderOnDiscreteFunctions() > 0)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 
      << "(*) discrete functions will require cofacet-based transformations");
    if (cellSpec==NoTermsNeedCofacets) 
    {
      cellSpec = SomeTermsNeedCofacets;
    }
  }
  
  SUNDANCE_MSG2(verb, tab << "found: " << cellSpec);
  
  return cellSpec;
}
  

void Assembler::configureVector(Array<Vector<double> >& b) const 
{
  /* Start timer, stopped upon dtor */
  TimeMonitor timer(configTimer());

  Tabs tab0;
  SUNDANCE_LEVEL1("vector config", tab0 << "in Assembler::configureVector()");

  /* Get the vector spaces for each block of equations */
  Array<VectorSpace<double> > vs(eqn_->numVarBlocks());
  for (unsigned int i=0; i<eqn_->numVarBlocks(); i++)
  {
    vs[i] = privateRowSpace_[i]->vecSpace();
  }
  VectorSpace<double> rowSpace;
  
  /* If we have more than one block, we make a Cartesian product space containing
   * the spaces for each block */
  if (eqn_->numVarBlocks() > 1)
  {
    rowSpace = TSFExtended::productSpace(vs);
  }
  else /* Otherwise we have a single, monolithic vector space */
  {
    rowSpace = vs[0];
  }

  /* Create each vector in the multivector */
  for (unsigned int i=0; i<b.size(); i++)
  {
    /* Create the vector. Recall that the vector space is a factory used to 
     * create a vector of specified size and distribution */
    b[i] = rowSpace.createMember();

    /* If the vector is blocked, configure the blocks */
    if (!partitionBCs_ && eqn_->numVarBlocks() > 1)
    {
      /* configure the blocks */
      Vector<double> vecBlock;
      for (unsigned int br=0; br<eqn_->numVarBlocks(); br++)
      {
        configureVectorBlock(br, vecBlock);
        b[i].setBlock(br, vecBlock);
      }
    }
    else  
    {
      /* nothing to do here except check that the vector is loadable */
      if (!partitionBCs_)
      {
        TSFExtended::LoadableVector<double>* lv 
          = dynamic_cast<TSFExtended::LoadableVector<double>* >(b[i].ptr().get());
        
        TEST_FOR_EXCEPTION(lv == 0, RuntimeError,
          "vector is not loadable in Assembler::configureVector()");
      }
      else
      {
      }
    }
  }
  numConfiguredColumns_ = b.size();
}

void Assembler::configureVectorBlock(int br, Vector<double>& b) const 
{
  Tabs tab0;
  SUNDANCE_LEVEL2("vector config", tab0 << "in Assembler::configureVectorBlock()");
  VectorSpace<double> vecSpace = privateRowSpace_[br]->vecSpace();

  b = vecSpace.createMember();
  
  if (!partitionBCs_)
  {
    TSFExtended::LoadableVector<double>* lv 
      = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());
    
    TEST_FOR_EXCEPTION(lv == 0, RuntimeError,
      "vector block is not loadable "
      "in Assembler::configureVectorBlock()");
  }
}


void Assembler::configureMatrix(LinearOperator<double>& A,
  Array<Vector<double> >& b) const
{
  TimeMonitor timer(configTimer());
  
  if (matNeedsConfiguration_)
  {
    Tabs tab0;
    SUNDANCE_LEVEL1("matrix config", tab0 << "in Assembler::configureMatrix()");
    int nRowBlocks = rowMap_.size();
    int nColBlocks = colMap_.size();
    Array<Array<int> > isNonzero = findNonzeroBlocks();

    if (nRowBlocks==1 && nColBlocks==1)
    {
      configureMatrixBlock(0,0,A);
    }
    else
    {
      A = makeBlockOperator(solnVecSpace(), rowVecSpace());
      for (int br=0; br<nRowBlocks; br++)
      {
        for (int bc=0; bc<nColBlocks; bc++)
        {
          if (isNonzero[br][bc])
          {
            LinearOperator<double> matBlock;
            configureMatrixBlock(br, bc, matBlock);
            A.setBlock(br, bc, matBlock);
          }
        }
      }
      A.endBlockFill();
    }
    matNeedsConfiguration_ = false;
  }
  else
  {
    Tabs tab0;
    SUNDANCE_LEVEL1("matrix config", 
      tab0 << "Assembler::configureMatrix() not needed, proceeding to configure vector");
  }
  configureVector(b);
}

void Assembler::configureMatrixBlock(int br, int bc,
  LinearOperator<double>& A) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());

  SUNDANCE_LEVEL1("matrix config", tab << "in Assembler::configureMatrixBlock()");
  
  SUNDANCE_LEVEL2("matrix config", tab << "Assembler: num rows = " << rowMap()[br]->numDOFs());
  
  SUNDANCE_LEVEL2("matrix config", tab << "Assembler: num cols = " << colMap()[bc]->numDOFs());
  
  VectorSpace<double> rowSpace = privateRowSpace_[br]->vecSpace();
  VectorSpace<double> colSpace = privateColSpace_[bc]->vecSpace();

  RefCountPtr<MatrixFactory<double> > matFactory ;

  if (partitionBCs_)
  {
    matFactory = rcp(new PartitionedMatrixFactory(colSpace, lowestCol_[bc],
        isBCCol_[bc], remoteBCCols_[bc], colVecType_[bc], 
        rowSpace, lowestRow_[br], isBCRow_[br], rowVecType_[br]));
  }
  else
  {
    matFactory = rowVecType_[br].createMatrixFactory(colSpace, rowSpace);
  }

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(matFactory.get());

  CollectivelyConfigurableMatrixFactory* ccmf 
    = dynamic_cast<CollectivelyConfigurableMatrixFactory*>(matFactory.get());

  TEST_FOR_EXCEPTION(ccmf==0 && icmf==0, RuntimeError,
    "Neither incremental nor collective matrix structuring "
    "appears to be available");


  /* If collective structuring is the user preference, or if incremental
   * structuring is not supported, do collective structuring */
  if (false /* (icmf==0 || !matrixEliminatesRepeatedCols()) && ccmf != 0 */)
  {
    Tabs tab1;
    SUNDANCE_LEVEL2("matrix config", tab1 << "Assembler: doing collective matrix structuring...");
    Array<int> graphData;
    Array<int> nnzPerRow;
    Array<int> rowPtrs;
      
    using Teuchos::createVector;

    getGraph(br, bc, graphData, rowPtrs, nnzPerRow);
    ccmf->configure(lowestRow_[br], createVector(rowPtrs), createVector(nnzPerRow), createVector(graphData));
  }
  else
  {
    Tabs tab1;
    SUNDANCE_LEVEL2("matrix config", tab1 << "Assembler: doing incremental matrix structuring...");
    incrementalGetGraph(br, bc, icmf);
    {
      TimeMonitor timer1(matFinalizeTimer());
      icmf->finalize();
    }
  }
  
  SUNDANCE_LEVEL3("matrix config", tab << "Assembler: allocating matrix...");
  {
    TimeMonitor timer1(matAllocTimer());
    A = matFactory->createMatrix();
  }
}

TSFExtended::LinearOperator<double> Assembler::allocateMatrix() const
{
  LinearOperator<double> A;
  Array<Vector<double> > b;
  configureMatrix(A, b);
  return A;
}





  
void Assembler::displayEvaluationResults(
  const EvalContext& context, 
  const EvaluatableExpr* evalExpr, 
  const Array<double>& constantCoeffs, 
  const Array<RefCountPtr<EvalVector> >& vectorCoeffs) const 
{
  Tabs tab;
  FancyOStream& os = Out::os();

  os << tab << "evaluation results: " << endl;

  const RefCountPtr<SparsitySuperset>& sparsity 
    = evalExpr->sparsitySuperset(context);
  
  sparsity->print(os, vectorCoeffs, constantCoeffs);
}



void Assembler::assemblyLoop(const ComputationType& compType,
  RefCountPtr<AssemblyKernelBase> kernel) const
{
  Tabs tab;
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);

  SUNDANCE_BANNER1(verb, tab, "Assembly loop");

  SUNDANCE_MSG2(verb, tab << "computation type is " << compType); 
  /* Allocate space for the workset's list of cell local IDs.
   * Currently, a workset is an array of cell indices. It could be an array
   * of pointers to cell objects, cell iterators, or whatever is needed to
   * work with something like a Peano curve data structure. */
  SUNDANCE_MSG2(verb, tab << "work set size is " << workSetSize()); 
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  /* Allocate isLocalFlag array, which distinguishes between local
   * and non-local cells in the workset. This is used to prevent 
   * adding multiple copies of zero-form values for border cells. */
  RefCountPtr<Array<int> > isLocalFlag = rcp(new Array<int>());
  isLocalFlag->reserve(workSetSize());

  /* Declare objects for storage of coefficient evaluation results 
   * that are returned from the symbolic evaluation system. EvalVector
   * is the object in which a vector of numerical values for a 
   * given functional derivative are returned from an evaluation of
   * a symbolic DAG. */
  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;

  /* Create an object in which to store local integration results */
  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  /* Get the symbolic specification of the current computation.
   * The "context" is simply a unique ID used to distinguish different
   * settings in which evaluation might be made. The same expression might be
   * evaluated in a context where some subset of all possible functional
   * derivatives are needed, and later in another context where a different 
   * subset of functional derivatives is needed. The "context" object lets
   * us associate a different Evaluator object with each such set of
   * requirements. */
  const Array<EvalContext>& contexts = contexts_.get(compType);
  /* Get the integral groups needed for this calculation */
  const Array<Array<RCP<IntegralGroup> > >& groups = groups_.get(compType);
  /* Get the expressions needed for this calculation */
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(compType);

  
  /* === Start main loop. 
   * The outer loop is over RQCs, that is, over unique combinations of subregions
   * (CellFilters) and quadrature rules.   */
  SUNDANCE_MSG1(verb, 
    tab << "---------- outer assembly loop over subregions");

  /* Record the default kernel verbosity so that it if changes we can
   * reset it at the end of a loop iteration */
  int oldKernelVerb = kernel->verb();
  
  /* Looping over RQCs */
  for (unsigned int r=0; r<rqc_.size(); r++)
  {
      Tabs tab0;

      /* Set the verbosity level for this RQC */
      int rqcVerb = verb;
      int evalVerb = 0;

      /* Check for watch point, and set verbosity accordingly */
      if (rqc_[r].watch().isActive()) 
      {
        rqcVerb=5;
        evalVerb = rqc_[r].watch().param("evaluation");
        kernel->setVerbosity(rqc_[r].watch().param("fill"));

        SUNDANCE_MSG2(rqcVerb, tab0 << endl 
          << tab0 << "-------------"
          << endl << tab0 << " doing watched subregion=" 
          << rqc_[r]);     
      }
      Tabs tab01;

      SUNDANCE_MSG2(rqcVerb, tab01 << "expr is " << evalExprs[r]->toString());
      SUNDANCE_MSG2(rqcVerb, tab01 << "isBC= " << isBCRqc_[r]);

    
      /* Deciding whether we should skip this RQC in the current computation 
       * type. For example, a certain boundary surface might appear in the
       * computation of a functional but not in the state equation. */
      if ((!isBCRqc_[r] && eqn_->skipRqc(compType, rqc_[r]))
        || (isBCRqc_[r] && eqn_->skipBCRqc(compType, rqc_[r])))
      {
        SUNDANCE_MSG2(rqcVerb, "nothing to do for comp type " 
          << compType);
        continue;
      }

      /* specify the evaluation mediator for this RQC.
       * Recall that the evaluation mediator is the object responsible for communication
       * between the symbolic expression tree and discretization-dependent data structures
       * such as discrete functions. 
       *
       * The evaluation manager is an object that is responsible for management
       * of the symbolic evaluation; it controls allocation of memory for evaluation
       * results, access to the evaluation mediator, and other administrative tasks. 
       */
      evalMgr_->setMediator(mediators_[r]);
      /* Tell the manager which CellFilter and QuadratureFamily we're currently working with. 
       * This is simply forwarded to the mediator, which needs to know the number
       * of quadrature points as well as mesh properties such as cell dimension. */
      evalMgr_->setRegion(contexts_.get(compType)[r]);
  
      /* get the cell filter for the current RQC */
      CellFilter filter = rqc_[r].domain();
      /* Find the cells that "pass" the predicate in the filter. Note: a CellFilter
       * will cache the cell sets it has computed, so the predicate computation will
       * only be done once per mesh, regardless of how often this function is called. 
       * Todo: the cache needs to be reset upon mesh refinement. 
       */
      CellSet cells = filter.getCells(mesh_);
      /* Find the dimension of cells that pass the current filter. */
      int cellDim = filter.dimension(mesh_);
      /* Find the type of cells in the current filter. Note: we've assumed here that all
       * cells have identical topology, and need to work out how to deal with 
       * meshes with mxed cell types. That will usually be handled by grouping similar
       * cells into "blocks" (as is done in Exodus files, for instance) in which case 
       * will still work, but there should be an error check to ensure that that assumption
       * is never violated. */ 
      CellType cellType = mesh_.cellType(cellDim);
      /* Get the cell type of the maximal-dimension cofacets, in case we need 
       * integrals or DOF maps done on the maximal cofacets. Todo: this code will break
       * for internal boundaries at the interface between cofacets of different types,
       * e.g., a face joining a prism and a hex. Not sure how to handle that case. */
      CellType maxCellType = mesh_.cellType(mesh_.spatialDim());

      SUNDANCE_MSG2(rqcVerb, tab01 << "cell type = " << cellType 
        << ", cellDim = " << cellDim 
        << ", max cell type = " << maxCellType 
        << ", max cell dim = " << mesh_.spatialDim());


      /* Determine whether we need to refer to maximal cofacets for 
       * some or all integrations and DOF mappings. */
      IntegrationCellSpecifier intCellSpec
        = rqcRequiresMaximalCofacets_.get(compType)[r];
      SUNDANCE_MSG2(rqcVerb, tab01 
        << "whether we need to refer to maximal cofacets: " << intCellSpec);

      /* Find the unknowns and variations appearing on the current domain. This
       * information is stored in the EquationSet object.  */
      const Array<Set<int> >& requiredVars = eqn_->reducedVarsOnRegion(filter);
      const Array<Set<int> >& requiredUnks = eqn_->reducedUnksOnRegion(filter);

      /* Prepare for evaluation on the current domain:
       * Tell the mediator whether maximal cofacets should be used (which will determine
       * how discrete functions are computed). */
      mediators_[r]->setIntegrationSpec(intCellSpec);
      /*
       * Tell the mediator the cell type, and whether we are on an internal
       * boundary. We need to know if we're on an internal boundary so that 
       * we can use DOF lookup logic that's capable of figuring out which
       * of two cofacets contains DOF information for those functions defined
       * on only one side of the boundary (as in, e.g., fluid-structure boundaries).
       */
      mediators_[r]->setCellType(cellType, maxCellType, isInternalBdry_[r]);    
      /* Get the Evaluator object that will actually carry out calculations on
       * the expression DAG in the current context (recall that a single expression
       * may support multiple evaluators). */
      const Evaluator* evaluator 
        = evalExprs[r]->evaluator(contexts[r]).get();

      /* Loop over cells in batches of the work set size.
       * At present, we're accumulating cell indices into an array. That would
       * need to be changed to work with Peano. */
      CellIterator iter=cells.begin();
      int workSetCounter = 0;
      int myRank = mesh_.comm().getRank();

      SUNDANCE_MSG2(rqcVerb, tab01 << "----- looping over worksets");
      while (iter != cells.end())
      {
        Tabs tab1;
        /* build up the work set: add cells until the work set size is 
         * reached or we run out of cells. To begin with, empty both the
         * workset array and the isLocalFlag array. Note that the reserve()
         * method has been called previously, so that as we append cells
         * to the array, no memory allocation is done (unless we run over 
         * the reserved size). */
        workSet->resize(0);
        isLocalFlag->resize(0);
        for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
        {
          workSet->append(*iter);
          /* we need the isLocalFlag values so that we can ignore contributions
           * to zero-forms from off-processor elements */
          isLocalFlag->append(myRank==mesh_.ownerProcID(cellDim, *iter));
        }
        /* The work set has now been accumulated */
        SUNDANCE_MSG2(rqcVerb,
          tab1 << "doing work set " << workSetCounter
          << " consisting of " << workSet->size() << " cells");
        {
          Tabs tab2;
          SUNDANCE_MSG4(rqcVerb, tab2 << "cells are " << *workSet);
        }
        workSetCounter++;

        /* set the verbosity for the evaluation mediator */
        if (evalVerb > 0)
        { 
          evalMgr_->setVerbosity(evalVerb);
        }

        /* Register the workset with the mediator. Internally, the mediator
         * will look up the cell Jacobians and facet indices needed for this calculation. It 
         * uses them for discrete function transformation. */
        mediators_[r]->setCellBatch(workSet);
        /* Get the Jacobians from the mediator, so that we can use the same Jacobian
         * objects for discrete function transformation and for integral 
         * transformation. The "volume" Jacobian is used to scale the integration
         * cell volume by det(J). The "transformation" Jacobian is used to
         * transform vectors. These will be the same, except for the case
         * where we integrate on a facet but do transformations by reference to 
         * a maximal cofacet. */
        const CellJacobianBatch& JVol = mediators_[r]->JVol();
        const CellJacobianBatch& JTrans = mediators_[r]->JTrans();
        /* Get from the mediator the facet indices for each cell. If I am integrating
         * on a facet (e.g., a boundary cell) but getting DOFs or JTrans from
         * a maximal cofacet, I need to know my index in the array of 
         * that cofacet's facets. */
        const Array<int>& facetIndices = mediators_[r]->facetIndices();

        /* Reset the assembly kernel for the current workset. What happens at this
         * step depends on the specific kernel being used. The kernel might, for instance,
         * build local DOF maps for the current batch of cells. */
        kernel->prepareForWorkSet(
          requiredVars, 
          requiredUnks, 
          mediators_[r]);
        
        /* Evaluate the coefficient expressions. Recall that each coefficient
         * appearing in an integral is a particular functional derivative of the
         * integrand. The constant-valued coefficients are written into the
         * constantCoeffs array, the variable coefficients into the vectorCoeffs
         * array. 
         *
         * Recall that the eval manager contains the current evaluation mediator, so that
         * by passing the evaluation manager as an argument to evaluate(), the evaluation
         * is aware of the mediator and can therefore access discrete functions, etc.
         */ 
        evaluator->resetNumCalls();
        SUNDANCE_MSG2(rqcVerb, tab1 
          << "====== evaluating coefficient expressions") ;
        try
        {
          evalExprs[r]->evaluate(*evalMgr_, constantCoeffs, vectorCoeffs);
        }
        catch(std::exception& exc)
        {
          Tabs tabX;
          SUNDANCE_BANNER1(10, tabX, 
            "DETECTED EXCEPTION DURING EXPR EVALUATION CALL IN ASSEMBLY LOOP");
          Tabs tabX1;
          SUNDANCE_MSG1(10, tabX1 << "While working on RQC="
            << rqc_[r]);
          SUNDANCE_MSG1(10, tabX1 << "While evaluating expr="
            << evalExprs[r]->toString());
          throw (RuntimeError(exc.what()));
        }

        /* Optionally, print the evaluation results */
        SUNDANCE_MSG2(rqcVerb, tab1 << "====== done evaluating expressions") ;
        if (evalVerb > 2) 
        {
          displayEvaluationResults(contexts[r], evalExprs[r], constantCoeffs, 
            vectorCoeffs);
        }
    
        /* ---- Do the element integrals and insertions ------ */

        /* The matrices used to transform integrals are built upon first use by 
         * this workset, then cached because they may be needed for several 
         * integrals on the same workset. As we're now in a new workset with
         * new cells, they should be rebuilt if needed. This step informs all integrals
         * that the cache is out of date. 
         *
         * Todo: this uses a static function that contains static data (via the "Meyers
         * trick") and so is not thread-safe. If we want to do multithreaded parallelism
         * for multicore architectures, this implementation will need to be changed.
         */ 
        ElementIntegral::invalidateTransformationMatrices();
      
        /* Loop over the integral groups */
        SUNDANCE_MSG2(rqcVerb, tab1 << "----- looping over integral groups");
        for (unsigned int g=0; g<groups[r].size(); g++)
        {
          Tabs tab2;
          SUNDANCE_MSG2(rqcVerb, tab2 << endl << tab2 
            << "--- evaluating integral group g=" << g << " of " 
            << groups[r].size() );

          /* Do the integrals. The integration results will be written into
           * the array "localValues". */
          const RCP<IntegralGroup>& group = groups[r][g];
          if (!group->evaluate(JTrans, JVol, *isLocalFlag, facetIndices, 
              vectorCoeffs, constantCoeffs, localValues)) continue;
          /* add the integration results into the output objects by a call
           * to the kernel's fill() function. We need to pass isBCRqc to the kernel
           * because it might handle BC rows differently. The integral group
           * data structure contains information about which test and unknown
           * functions were used in this integration, and so provides to the assembly
           * kernel such information as is needed to look up the correct DOFs for this
           * batch of integrals. */
          kernel->fill(isBCRqc_[r], *group, localValues);
        }
        SUNDANCE_MSG2(rqcVerb, tab1 << "----- done looping over integral groups");
      }
      SUNDANCE_MSG2(rqcVerb, tab0 << "----- done looping over worksets");
      /* reset the kernel verbosity to the default */
      kernel->setVerbosity(oldKernelVerb);
  }
  SUNDANCE_MSG2(verb, tab << "----- done looping over rqcs");


  /* Do any post-fill processing, such as MPI_AllReduce add on functional values. */
  SUNDANCE_MSG2(verb, tab << "doing post-loop processing"); 
  kernel->postLoopFinalization();
  SUNDANCE_BANNER1(verb, tab, "done assembly loop"); 

  /* All done with assembly! */
}



/* ------------  assemble both the vector and the matrix  ------------- */

void Assembler::assemble(LinearOperator<double>& A,
  Array<Vector<double> >& mv) const 
{
  TimeMonitor timer(assemblyTimer());
  Tabs tab;
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);
  
  SUNDANCE_BANNER1(verb, tab, "Assembling matrix and vector");

  TEST_FOR_EXCEPTION(!contexts_.containsKey(MatrixAndVector),
    RuntimeError,
    "Assembler::assemble(A, b) called for an assembler that "
    "does not support matrix/vector assembly");

  configureMatrix(A, mv);
  
  RefCountPtr<AssemblyKernelBase> kernel 
    = rcp(new MatrixVectorAssemblyKernel(
            rowMap_, isBCRow_, lowestRow_,
            colMap_, isBCCol_, lowestCol_,
            A, mv, partitionBCs_, 
            0));

  assemblyLoop(MatrixAndVector, kernel);

  SUNDANCE_MSG1(verb, tab << "matrix=" << A);
  if (verb>0) A.print(Out::os());
  SUNDANCE_MSG1(verb, tab << "vectors=" << mv);
  for (unsigned int i=0; i<mv.size(); i++) 
  {
    SUNDANCE_MSG1(verb, tab << "vectors #" << i);
    if (verb>0) mv[i].print(Out::os());
  }

  SUNDANCE_MSG1(verb, tab << "Assembler: done assembling matrix and vector");
}

/* ------------  assemble the matrix and sensitivity RHS ------------- */

void Assembler::assembleSensitivities(LinearOperator<double>& A,
  Array<Vector<double> >& mv) const 
{
  TimeMonitor timer(assemblyTimer());
  Tabs tab;
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);
  
  SUNDANCE_BANNER1(verb, tab, "Assembling matrix and sensitivity vector");

  TEST_FOR_EXCEPTION(!contexts_.containsKey(Sensitivities),
    RuntimeError,
    "Assembler::assembleSensitivities(A, b) called for an assembler that "
    "does not support sensitivity assembly");

  configureMatrix(A, mv);
  
  
  RefCountPtr<AssemblyKernelBase> kernel 
    = rcp(new MatrixVectorAssemblyKernel(
            rowMap_, isBCRow_, lowestRow_,
            colMap_, isBCCol_, lowestCol_,
            A, mv, partitionBCs_, 
            0));

  assemblyLoop(Sensitivities, kernel);
  SUNDANCE_MSG1(verb, tab << "Assembler: done assembling matrix and sensitivity vector");
}


/* ------------  assemble the vector alone  ------------- */

void Assembler::assemble(Array<Vector<double> >& mv) const 
{
  /* Tab is advanced by ctor, taken back by dtor upon leaving scope */
  Tabs tab;
  /* Timer is started by ctor, stopped by dtor upon leaving scope */
  TimeMonitor timer(assemblyTimer());  

  /* If any subexpression is watched, print basic information */ 
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);

  SUNDANCE_BANNER1(verb, tab, "Assembling vector");

  /* Throw an exception if we don't know how to compute a vector */
  TEST_FOR_EXCEPTION(!contexts_.containsKey(VectorOnly),
    RuntimeError,
    "Assembler::assemble(b) called for an assembler that "
    "does not support vector-only assembly");

  /* The vector is not configured yet. Do so here. */
  configureVector(mv);
  
  /* Create an assembly kernel that knows how to fill a vector */
  RefCountPtr<AssemblyKernelBase> kernel 
    = rcp(new VectorAssemblyKernel(
            rowMap_, isBCRow_, lowestRow_,
            mv, partitionBCs_, 0));

  assemblyLoop(VectorOnly, kernel);

  SUNDANCE_MSG1(verb, tab << "Assembler: done assembling vector");
}


/* ------------  evaluate a functional and its gradient ---- */

void Assembler::evaluate(double& value, Array<Vector<double> >& gradient) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);

  SUNDANCE_BANNER1(verb, tab, "Computing functional and gradient");

  TEST_FOR_EXCEPTION(!contexts_.containsKey(FunctionalAndGradient),
    RuntimeError,
    "Assembler::evaluate(f,df) called for an assembler that "
    "does not support value/gradient assembly");

  configureVector(gradient);

  value = 0.0;
  
  RefCountPtr<AssemblyKernelBase> kernel 
    = rcp(new FunctionalGradientAssemblyKernel(
            mesh_.comm(),
            rowMap_, isBCRow_, lowestRow_,
            gradient, partitionBCs_, &value, verb));

  assemblyLoop(FunctionalAndGradient, kernel);

  SUNDANCE_BANNER1(verb, tab, "Done computing functional and gradient");

}




/* ------------  evaluate a functional ---- */

void Assembler::evaluate(double& value) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  int verb = 0;
  if (eqn_->hasActiveWatchFlag()) verb = max(verb, 1);

  SUNDANCE_BANNER1(verb, tab, "Computing functional");

  TEST_FOR_EXCEPTION(!contexts_.containsKey(FunctionalOnly),
    RuntimeError,
    "Assembler::evaluate(f) called for an assembler that "
    "does not support functional evaluation");

  value = 0.0;
  
  RefCountPtr<AssemblyKernelBase> kernel 
    = rcp(new FunctionalAssemblyKernel(mesh_.comm(), &value, verb));

  assemblyLoop(FunctionalOnly, kernel);

  SUNDANCE_BANNER1(verb, tab, "Done computing functional");
}




/* ------------  get the nonzero pattern for the matrix ------------- */
                       
                       
void Assembler::getGraph(int br, int bc, 
  Array<int>& graphData,
  Array<int>& rowPtrs,
  Array<int>& nnzPerRow) const 
{
  TimeMonitor timer(graphBuildTimer());
  Tabs tab;




  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  RefCountPtr<Array<Array<int> > > testLocalDOFs 
    = rcp(new Array<Array<int> >());

  RefCountPtr<Array<Array<int> > > unkLocalDOFs
    = rcp(new Array<Array<int> >());

  SUNDANCE_VERB_HIGH(tab << "Creating graph: there are " << rowMap_[br]->numLocalDOFs()
    << " local equations");


  Array<Set<int> > tmpGraph;
  tmpGraph.resize(rowMap_[br]->numLocalDOFs());

  {
    TimeMonitor timer2(colSearchTimer());
    for (unsigned int d=0; d<eqn_->numRegions(); d++)
    {
      Tabs tab0;
      CellFilter domain = eqn_->region(d);
      const Array<Set<int> >& requiredVars = eqn_->reducedVarsOnRegion(domain);
      const Array<Set<int> >& requiredUnks = eqn_->reducedUnksOnRegion(domain);
      SUNDANCE_OUT(this->verb() > 2, 
        tab0 << "cell set " << domain
        << " isBCRegion=" << eqn_->isBCRegion(d));
      unsigned int dim = domain.dimension(mesh_);
      CellSet cells = domain.getCells(mesh_);

      RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
      if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

      SUNDANCE_OUT(this->verb() > 2 && pairs.get() != 0, 
        tab0 << "non-BC pairs = "
        << *pairs);
       
      RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
      if (eqn_->isBCRegion(d))
      {
        if (eqn_->hasBCVarUnkPairs(domain)) 
        {
          bcPairs = eqn_->bcVarUnkPairs(domain);
          SUNDANCE_OUT(this->verb() > 2, tab0 << "BC pairs = "
            << *bcPairs);
        }
      }
      Array<Set<int> > unksForTestsSet(eqn_->numVarIDs(bc));
      Array<Set<int> > bcUnksForTestsSet(eqn_->numVarIDs(bc));

      Set<OrderedPair<int, int> >::const_iterator i;
      
      if (pairs.get() != 0)
      {
        for (i=pairs->begin(); i!=pairs->end(); i++)
        {
          const OrderedPair<int, int>& p = *i;
          int t = p.first();
          int u = p.second();

          TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
            "Test function ID " << t << " does not appear "
            "in equation set");
          TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
            "Unk function ID " << u << " does not appear "
            "in equation set");

          if (eqn_->blockForVarID(t) != br) continue;
          if (eqn_->blockForUnkID(u) != bc) continue;

          unksForTestsSet[eqn_->reducedVarID(t)].put(eqn_->reducedUnkID(u));
        }
      }
      if (bcPairs.get() != 0)
      {
        for (i=bcPairs->begin(); i!=bcPairs->end(); i++)
        {
          const OrderedPair<int, int>& p = *i;
          int t = p.first();
          int u = p.second();

          if (eqn_->blockForVarID(t) != br) continue;
          if (eqn_->blockForUnkID(u) != bc) continue;

          TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
            "Test function ID " << t << " does not appear "
            "in equation set");
          TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
            "Unk function ID " << u << " does not appear "
            "in equation set");
          bcUnksForTestsSet[eqn_->reducedVarID(t)].put(eqn_->reducedUnkID(u));
        }
      }

      Array<Array<int> > unksForTests(unksForTestsSet.size());
      Array<Array<int> > bcUnksForTests(bcUnksForTestsSet.size());

      for (unsigned int t=0; t<unksForTests.size(); t++)
      {
        unksForTests[t] = unksForTestsSet[t].elements();
        bcUnksForTests[t] = bcUnksForTestsSet[t].elements();
      }
      
      Array<int> numTestNodes;
      Array<int> numUnkNodes;

      

      int highestRow = lowestRow_[br] + rowMap_[br]->numLocalDOFs();

      int nt = eqn_->numVarIDs(br);

      CellIterator iter=cells.begin();
      while (iter != cells.end())
      {
        /* build a work set */
        workSet->resize(0);
        for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
        {
          workSet->append(*iter);
        }

        int nCells = workSet->size();

        RefCountPtr<const MapStructure> colMapStruct; 
        RefCountPtr<const MapStructure> rowMapStruct 
          = rowMap_[br]->getDOFsForCellBatch(dim, *workSet, 
            requiredVars[br], *testLocalDOFs,
            numTestNodes, this->verb());
        if (rowMap_[br].get()==colMap_[bc].get())
        {
          unkLocalDOFs = testLocalDOFs;
          numUnkNodes = numTestNodes;
          colMapStruct = rowMapStruct;
        }
        else
        {
          colMapStruct = colMap_[br]->getDOFsForCellBatch(dim, *workSet, 
            requiredUnks[bc], 
            *unkLocalDOFs, numUnkNodes, this->verb());
        }

        if (pairs.get() != 0)
        {
          for (int c=0; c<nCells; c++)
          {
            for (int t=0; t<nt; t++)
            {
              int tChunk = rowMapStruct->chunkForFuncID(t);
              int nTestFuncs = rowMapStruct->numFuncs(tChunk);
              int testFuncIndex = rowMapStruct->indexForFuncID(t);
              int nTestNodes = numTestNodes[tChunk];
              const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
              for (unsigned int uit=0; uit<unksForTests[t].size(); uit++)
              {
                Tabs tab2;
                int u = unksForTests[t][uit];
                int uChunk = colMapStruct->chunkForFuncID(u);
                int nUnkFuncs = colMapStruct->numFuncs(uChunk);
                int unkFuncIndex = colMapStruct->indexForFuncID(u);
                const Array<int>& unkDOFs = (*unkLocalDOFs)[uChunk];
                int nUnkNodes = numUnkNodes[uChunk];
                for (int n=0; n<nTestNodes; n++)
                {
                  int row
                    = testDOFs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
                  if (row < lowestRow_[br] || row >= highestRow
                    || (*(isBCRow_[br]))[row-lowestRow_[br]]) continue;
                  Set<int>& colSet = tmpGraph[row-lowestRow_[br]];
                  for (int m=0; m<nUnkNodes; m++)
                  {
                    int col 
                      = unkDOFs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes + m];
                    colSet.put(col);
                  }
                }
              }
            }
          }
        }
        if (bcPairs.get() != 0)
        {
          for (int c=0; c<nCells; c++)
          {
            for (int t=0; t<nt; t++)
            {
              int tChunk = rowMapStruct->chunkForFuncID(t);
              int nTestFuncs = rowMapStruct->numFuncs(tChunk);
              int testFuncIndex = rowMapStruct->indexForFuncID(t);
              int nTestNodes = numTestNodes[tChunk];
              const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
              for (unsigned int uit=0; uit<bcUnksForTests[t].size(); uit++)
              {
                Tabs tab2;
                int u = bcUnksForTests[t][uit];
                int uChunk = colMapStruct->chunkForFuncID(u);
                int nUnkFuncs = colMapStruct->numFuncs(uChunk);
                int unkFuncIndex = colMapStruct->indexForFuncID(u);
                const Array<int>& unkDOFs = (*unkLocalDOFs)[uChunk];
                int nUnkNodes = numUnkNodes[uChunk];
                for (int n=0; n<nTestNodes; n++)
                {
                  int row
                    = testDOFs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
                  if (row < lowestRow_[br] || row >= highestRow
                    || !(*(isBCRow_[br]))[row-lowestRow_[br]]) continue;
                  Set<int>& colSet = tmpGraph[row-lowestRow_[br]];
                  for (int m=0; m<nUnkNodes; m++)
                  {
                    int col 
                      = unkDOFs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes + m];
                    colSet.put(col);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  
  {
    TimeMonitor t2(graphFlatteningTimer());
    unsigned int nLocalRows = rowMap_[br]->numLocalDOFs();

    unsigned int nnz = 0;
    rowPtrs.resize(nLocalRows);
    nnzPerRow.resize(rowMap_[br]->numLocalDOFs());
    for (unsigned int i=0; i<nLocalRows; i++) 
    {
      rowPtrs[i] = nnz;
      nnzPerRow[i] = tmpGraph[i].size();
      nnz += nnzPerRow[i];
    }

    graphData.resize(nnz);
    int* base = &(graphData[0]);
    for (unsigned int i=0; i<nLocalRows; i++)
    {
      //        tmpGraph[i].fillArray(base + rowPtrs[i]);
      int* rowBase = base + rowPtrs[i];
      const Set<int>& rowSet = tmpGraph[i];
      int k = 0;
      for (Set<int>::const_iterator 
             j=rowSet.begin(); j != rowSet.end(); j++, k++)
      {
        rowBase[k] = *j;
      }
    }
  }

}
/* ------------  get the nonzero pattern for the matrix ------------- */
                       
                       
void Assembler
::incrementalGetGraph(int br, int bc,
  IncrementallyConfigurableMatrixFactory* icmf) const 
{
  TimeMonitor timer(graphBuildTimer());
  Tabs tab;


  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  RefCountPtr<Array<Array<int> > > testLocalDOFs 
    = rcp(new Array<Array<int> >());

  RefCountPtr<Array<Array<int> > > unkLocalDOFs
    = rcp(new Array<Array<int> >());

  SUNDANCE_VERB_LOW(tab << "Creating graph: there are " << rowMap_[br]->numLocalDOFs()
    << " local equations");


  for (unsigned int d=0; d<eqn_->numRegions(); d++)
  {
    Tabs tab0;
    CellFilter domain = eqn_->region(d);
    const Array<Set<int> >& requiredVars = eqn_->reducedVarsOnRegion(domain);
    const Array<Set<int> >& requiredUnks = eqn_->reducedUnksOnRegion(domain);
    Array<int> localVars = requiredVars[br].elements();
    Array<int> localUnks = requiredUnks[bc].elements();
    SUNDANCE_OUT(this->verb() > 2, 
      tab0 << "cell set " << domain
      << " isBCRegion=" << eqn_->isBCRegion(d));
    unsigned int dim = domain.dimension(mesh_);
    CellSet cells = domain.getCells(mesh_);

    RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
    if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

    SUNDANCE_OUT(this->verb() > 2 && pairs.get() != 0, 
      tab0 << "non-BC pairs = "
      << *pairs);
       
    RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
    if (eqn_->isBCRegion(d))
    {
      if (eqn_->hasBCVarUnkPairs(domain)) 
      {
        bcPairs = eqn_->bcVarUnkPairs(domain);
        SUNDANCE_OUT(this->verb() > 2, tab0 << "BC pairs = "
          << *bcPairs);
      }
    }
    Array<Set<int> > unksForTestsSet(eqn_->numVarIDs(br));
    Array<Set<int> > bcUnksForTestsSet(eqn_->numVarIDs(br));

    Set<OrderedPair<int, int> >::const_iterator i;
      
    if (pairs.get() != 0)
    {
      for (i=pairs->begin(); i!=pairs->end(); i++)
      {
        const OrderedPair<int, int>& p = *i;
        int t = p.first();
        int u = p.second();

        TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
          "Test function ID " << t << " does not appear "
          "in equation set");
        TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
          "Unk function ID " << u << " does not appear "
          "in equation set");


        if (eqn_->blockForVarID(t) != br) continue;
        if (eqn_->blockForUnkID(u) != bc) continue;

        unksForTestsSet[eqn_->reducedVarID(t)].put(eqn_->reducedUnkID(u));
      }
    }
    if (bcPairs.get() != 0)
    {
      for (i=bcPairs->begin(); i!=bcPairs->end(); i++)
      {
        const OrderedPair<int, int>& p = *i;
        int t = p.first();
        int u = p.second();
        TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
          "Test function ID " << t << " does not appear "
          "in equation set");
        TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
          "Unk function ID " << u << " does not appear "
          "in equation set");

        if (eqn_->blockForVarID(t) != br) continue;
        if (eqn_->blockForUnkID(u) != bc) continue;

        bcUnksForTestsSet[eqn_->reducedVarID(t)].put(eqn_->reducedUnkID(u));
      }
    }

    Array<Array<int> > unksForTests(unksForTestsSet.size());
    Array<Array<int> > bcUnksForTests(bcUnksForTestsSet.size());

    for (unsigned int t=0; t<unksForTests.size(); t++)
    {
      unksForTests[t] = unksForTestsSet[t].elements();
      bcUnksForTests[t] = bcUnksForTestsSet[t].elements();
    }
      
    Array<int> numTestNodes;
    Array<int> numUnkNodes;
      
    int highestRow = lowestRow_[br] + rowMap_[br]->numLocalDOFs();

    CellIterator iter=cells.begin();

    while (iter != cells.end())
    {
      /* build a work set */
      workSet->resize(0);
      for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
      {
        workSet->append(*iter);
      }

      int nCells = workSet->size();

      RefCountPtr<const MapStructure> colMapStruct; 
      RefCountPtr<const MapStructure> rowMapStruct 
        = rowMap_[br]->getDOFsForCellBatch(dim, *workSet, 
          requiredVars[br],
          *testLocalDOFs,
          numTestNodes, this->verb());

      if (rowMap_[br].get()==colMap_[bc].get())
      {
        unkLocalDOFs = testLocalDOFs;
        numUnkNodes = numTestNodes;
        colMapStruct = rowMapStruct;
      }
      else
      {
        colMapStruct = colMap_[bc]->getDOFsForCellBatch(dim, *workSet, 
          requiredUnks[bc],
          *unkLocalDOFs, numUnkNodes, this->verb());
      }

          
      if (pairs.get() != 0)
      {
        for (int c=0; c<nCells; c++)
        {
          for (unsigned int tIndex=0; tIndex<localVars.size(); tIndex++)
          {
            int t = localVars[tIndex];
            int tChunk = rowMapStruct->chunkForFuncID(t);
            int nTestFuncs = rowMapStruct->numFuncs(tChunk);
            int testFuncIndex = rowMapStruct->indexForFuncID(t);
            int nTestNodes = numTestNodes[tChunk];
            const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
            for (unsigned int uit=0; uit<unksForTests[t].size(); uit++)
            {
              Tabs tab2;
              int u = unksForTests[t][uit];
              int uChunk = colMapStruct->chunkForFuncID(u);
              int nUnkFuncs = colMapStruct->numFuncs(uChunk);
              int unkFuncIndex = colMapStruct->indexForFuncID(u);
              const Array<int>& unkDOFs = (*unkLocalDOFs)[uChunk];
              int nUnkNodes = numUnkNodes[uChunk];
              for (int n=0; n<nTestNodes; n++)
              {
                int row
                  = testDOFs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
                if (row < lowestRow_[br] || row >= highestRow
                  || (*(isBCRow_[br]))[row-lowestRow_[br]]) continue;
                const int* colPtr = &(unkDOFs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes]);
                icmf->initializeNonzerosInRow(row, nUnkNodes, colPtr);
              }
            }
          }
        }
      }
      if (bcPairs.get() != 0)
      {
        for (int c=0; c<nCells; c++)
        {
          for (unsigned int tIndex=0; tIndex<localVars.size(); tIndex++)
          {
            int t = localVars[tIndex];
            int tChunk = rowMapStruct->chunkForFuncID(t);
            int nTestFuncs = rowMapStruct->numFuncs(tChunk);
            int testFuncIndex = rowMapStruct->indexForFuncID(t);
            int nTestNodes = numTestNodes[tChunk];
            const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
            for (unsigned int uit=0; uit<bcUnksForTests[t].size(); uit++)
            {
              Tabs tab2;
              int u = bcUnksForTests[t][uit];
              int uChunk = colMapStruct->chunkForFuncID(u);
              int nUnkFuncs = colMapStruct->numFuncs(uChunk);
              int unkFuncIndex = colMapStruct->indexForFuncID(u);
              const Array<int>& unkDOFs = (*unkLocalDOFs)[uChunk];
              int nUnkNodes = numUnkNodes[uChunk];
              for (int n=0; n<nTestNodes; n++)
              {
                int row
                  = testDOFs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
                if (row < lowestRow_[br] || row >= highestRow
                  || !(*(isBCRow_[br]))[row-lowestRow_[br]]) continue;

                const int* colPtr = &(unkDOFs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes]);
                icmf->initializeNonzerosInRow(row, nUnkNodes, colPtr);
              }
            }
          }
        }
      }
    }
  }
}


Array<Array<int> > Assembler::findNonzeroBlocks() const
{
  Array<Array<int> > rtn(eqn_->numVarBlocks());
  for (unsigned int br=0; br<rtn.size(); br++)
  {
    rtn[br].resize(eqn_->numUnkBlocks());
    for (unsigned int bc=0; bc<rtn[br].size(); bc++)
    {
      rtn[br][bc] = 0 ;
    }
  }

  for (unsigned int d=0; d<eqn_->numRegions(); d++)
  {
    Tabs tab0;
    CellFilter domain = eqn_->region(d);
    SUNDANCE_OUT(this->verb() > 2, 
      tab0 << "cell set " << domain
      << " isBCRegion=" << eqn_->isBCRegion(d));

    RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
    if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

    SUNDANCE_OUT(this->verb() > 2 && pairs.get() != 0, 
      tab0 << "non-BC pairs = "
      << *pairs);
       
    RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
    if (eqn_->isBCRegion(d))
    {
      if (eqn_->hasBCVarUnkPairs(domain)) 
      {
        bcPairs = eqn_->bcVarUnkPairs(domain);
        SUNDANCE_OUT(this->verb() > 2, tab0 << "BC pairs = "
          << *bcPairs);
      }
    }

    Set<OrderedPair<int, int> >::const_iterator i;
      
    if (pairs.get() != 0)
    {
      for (i=pairs->begin(); i!=pairs->end(); i++)
      {
        const OrderedPair<int, int>& p = *i;
        int t = p.first();
        int u = p.second();

        TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
          "Test function ID " << t << " does not appear "
          "in equation set");
        TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
          "Unk function ID " << u << " does not appear "
          "in equation set");


        int br = eqn_->blockForVarID(t);
        int bc = eqn_->blockForUnkID(u);
        rtn[br][bc] = 1;
      }
    }
    if (bcPairs.get() != 0)
    {
      for (i=bcPairs->begin(); i!=bcPairs->end(); i++)
      {
        const OrderedPair<int, int>& p = *i;
        int t = p.first();
        int u = p.second();
        TEST_FOR_EXCEPTION(!eqn_->hasVarID(t), InternalError,
          "Test function ID " << t << " does not appear "
          "in equation set");
        TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
          "Unk function ID " << u << " does not appear "
          "in equation set");
        int br = eqn_->blockForVarID(t);
        int bc = eqn_->blockForUnkID(u);
        rtn[br][bc] = 1;
      }
    }
  }

  return rtn;
}

VectorSpace<double> Assembler::solnVecSpace() const
{
  Array<VectorSpace<double> > rtn(eqn_->numUnkBlocks());

  for (unsigned int i=0; i<rtn.size(); i++)
  {
    rtn[i] = solutionSpace()[i]->vecSpace();
  }

  if ((int) rtn.size() == 1)
  {
    return rtn[0];
  }
  return productSpace(rtn);
}


VectorSpace<double> Assembler::rowVecSpace() const
{
  Array<VectorSpace<double> > rtn(eqn_->numVarBlocks());

  for (unsigned int i=0; i<rtn.size(); i++)
  {
    rtn[i] = rowSpace()[i]->vecSpace();
  }

  if ((int) rtn.size() == 1)
  {
    return rtn[0];
  }
  return productSpace(rtn);
}


Vector<double> Assembler
::convertToMonolithicVector(const Array<Vector<double> >& internalBlock,
  const Array<Vector<double> >& bcBlock) const
{

  Array<VectorSpace<double> > spaces(bcBlock.size());
  Array<Vector<double> > v(bcBlock.size());

  SUNDANCE_CHECK_ARRAY_SIZE_MATCH(internalBlock, bcBlock);
  SUNDANCE_CHECK_ARRAY_SIZE_MATCH(internalBlock, privateColSpace_);

  for (unsigned int i=0; i<internalBlock.size(); i++)
  {
    VectorSpace<double> partSpace = privateColSpace_[i]->vecSpace();
    Vector<double> in = partSpace.createMember();
    in.setBlock(0, internalBlock[i]);
    in.setBlock(1, bcBlock[i]);
    Vector<double> out = externalColSpace_[i]->vecSpace().createMember();
    spaces[i] = externalColSpace_[i]->vecSpace();
    converter_[i]->convert(in, out);
    v[i] = out;
  }

  if (spaces.size() > 1) 
  {
    VectorSpace<double> rtnSpace = productSpace(spaces);
    Vector<double> rtn = rtnSpace.createMember();
    for (unsigned int i=0; i<spaces.size(); i++)
    {
      rtn.setBlock(i, v[i]);
    }
    return rtn;
  }
  else
  {
    return v[0];
  }
  
}



unsigned int& Assembler::workSetSize()
{
  static unsigned int rtn = defaultWorkSetSize(); 
  return rtn;
}
