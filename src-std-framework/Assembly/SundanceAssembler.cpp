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
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_HashTable.h"
#include "SundanceIntHashSet.hpp"
#include "TSFProductVectorSpace.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
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

static Time& matInsertTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("matrix insertion"); 
  return *rtn;
}

static Time& vecInsertTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("vector insertion"); 
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

static Time& tmpGraphBuildTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("tmp graph creation"); 
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


Assembler
::Assembler(const Mesh& mesh, 
            const RefCountPtr<EquationSet>& eqn,
            const Array<VectorType<double> >& rowVectorType,
            const Array<VectorType<double> >& colVectorType,
            const VerbositySetting& verb)
  : matNeedsConfiguration_(true),
    matNeedsFinalization_(true),
    vecNeedsConfiguration_(true),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    rowSpace_(eqn->numVarBlocks()),
    colSpace_(eqn->numUnkBlocks()),
    bcRows_(eqn->numVarBlocks()),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(eqn->numVarBlocks()),
    lowestRow_(eqn->numVarBlocks()),
    rowVecType_(rowVectorType),
    colVecType_(colVectorType),
    testIDToBlockMap_(),
    unkIDToBlockMap_()
{
  TimeMonitor timer(assemblerCtorTimer());
  verbosity() = verb;
  init(mesh, eqn);
}

Assembler
::Assembler(const Mesh& mesh, 
            const RefCountPtr<EquationSet>& eqn,
            const VerbositySetting& verb)
  : matNeedsConfiguration_(true),
    matNeedsFinalization_(true),
    vecNeedsConfiguration_(true),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    rowSpace_(eqn->numVarBlocks()),
    colSpace_(eqn->numUnkBlocks()),
    bcRows_(eqn->numVarBlocks()),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(eqn->numVarBlocks()),
    lowestRow_(eqn->numVarBlocks()),
    rowVecType_(),
    colVecType_(),
    testIDToBlockMap_(),
    unkIDToBlockMap_()
{
  TimeMonitor timer(assemblerCtorTimer());
  verbosity() = verb;
  init(mesh, eqn);
}

void Assembler::init(const Mesh& mesh, 
                     const RefCountPtr<EquationSet>& eqn)
{
  RefCountPtr<GrouperBase> grouper = rcp(new TrivialGrouper());
  grouper->verbosity() = verbosity();

  const Set<ComputationType>& compTypes = eqn->computationTypes();

  DOFMapBuilder mapBuilder;

  if (compTypes.contains(VectorOnly) 
      || compTypes.contains(FunctionalAndGradient))
    {
      mapBuilder = DOFMapBuilder(mesh, eqn);

      rowMap_ = mapBuilder.rowMap();
      isBCRow_ = mapBuilder.isBCRow();
      lowestRow_.resize(eqn_->numVarBlocks());
      for (unsigned int b=0; b<lowestRow_.size(); b++) 
        {
          lowestRow_[b] = rowMap_[b]->lowestLocalDOF();
          rowSpace_[b] = rcp(new DiscreteSpace(mesh, mapBuilder.testBasisArray()[b], 
                                               rowMap_[b], rowVecType_[b]));
        }
    }

  if (!eqn->isFunctionalCalculator())
    {
      colMap_ = mapBuilder.colMap();
      for (unsigned int b=0; b<eqn_->numUnkBlocks(); b++) 
        {
          colSpace_[b] = rcp(new DiscreteSpace(mesh, mapBuilder.unkBasisArray()[b], 
                                               colMap_[b], colVecType_[b]));
        }
      groups_.put(MatrixAndVector, Array<Array<IntegralGroup> >());
      rqcRequiresMaximalCofacets_.put(MatrixAndVector, Array<int>());
      contexts_.put(MatrixAndVector, Array<EvalContext>());
      evalExprs_.put(MatrixAndVector, Array<const EvaluatableExpr*>());
      groups_.put(VectorOnly, Array<Array<IntegralGroup> >());
      rqcRequiresMaximalCofacets_.put(VectorOnly, Array<int>());
      contexts_.put(VectorOnly, Array<EvalContext>());
      evalExprs_.put(VectorOnly, Array<const EvaluatableExpr*>());
    }
  else
    {
      groups_.put(FunctionalAndGradient, Array<Array<IntegralGroup> >());
      rqcRequiresMaximalCofacets_.put(FunctionalAndGradient, Array<int>());
      contexts_.put(FunctionalAndGradient, Array<EvalContext>());
      evalExprs_.put(FunctionalAndGradient, Array<const EvaluatableExpr*>());
      groups_.put(FunctionalOnly, Array<Array<IntegralGroup> >());
      rqcRequiresMaximalCofacets_.put(FunctionalOnly, Array<int>());
      contexts_.put(FunctionalOnly, Array<EvalContext>());
      evalExprs_.put(FunctionalOnly, Array<const EvaluatableExpr*>());
    }

  for (unsigned int r=0; r<eqn->regionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];
                         
      rqc_.append(rqc);
      isBCRqc_.append(false);
      const Expr& expr = eqn->expr(rqc);

      SUNDANCE_VERB_HIGH("creating integral groups for rqc=" << rqc << endl
                         << "expr = " << expr);

      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      CellType maxCellType = mesh.cellType(mesh.spatialDim());
      QuadratureFamily quad(rqc.quad());

      for (Set<ComputationType>::const_iterator 
             i=eqn->computationTypes().begin(); 
           i!=eqn->computationTypes().end();
           i++)
        {
          const ComputationType& compType = *i;
          //          const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(compType, rqc);
          EvalContext context = eqn->rqcToContext(compType, rqc);
          contexts_[compType].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[compType].append(ee);
          const RefCountPtr<SparsitySuperset>& sparsity 
            = ee->sparsitySuperset(context);
          SUNDANCE_VERB_EXTREME("sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, maxCellType, mesh.spatialDim(),
                              cellType, cellDim, quad, sparsity, groups);
          groups_[compType].append(groups);
          bool reqCofacets = false;
          for (unsigned int g=0; g<groups.size(); g++)
            {
              if (groups[g].requiresMaximalCofacet()) 
                {
                  reqCofacets = true;
                  break;
                }
            }
          rqcRequiresMaximalCofacets_[compType].append(reqCofacets);
        }
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }

  
  for (unsigned int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];
      rqc_.append(rqc);
      isBCRqc_.append(true);
      const Expr& expr = eqn->bcExpr(rqc);

      SUNDANCE_VERB_HIGH("creating integral groups for rqc=" << rqc << endl
                         << "expr = " << expr.toXML().toString());
      
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      CellType maxCellType = mesh.cellType(mesh.spatialDim());
      QuadratureFamily quad(rqc.quad());

      for (Set<ComputationType>::const_iterator 
             i=eqn->computationTypes().begin(); 
           i!=eqn->computationTypes().end();
           i++)
        {
          const ComputationType& compType = *i;
          //        const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(compType, rqc);
          EvalContext context = eqn->bcRqcToContext(compType, rqc);
          contexts_[compType].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[compType].append(ee);
          const RefCountPtr<SparsitySuperset>& sparsity 
            = ee->sparsitySuperset(context);
          SUNDANCE_VERB_EXTREME("sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, maxCellType, mesh.spatialDim(),
                              cellType, cellDim, quad, sparsity, groups);
          groups_[compType].append(groups);
          bool reqCofacets = false;
          for (unsigned int g=0; g<groups.size(); g++)
            {
              if (groups[g].requiresMaximalCofacet()) 
                {
                  reqCofacets = true;
                  break;
                }
            }
          rqcRequiresMaximalCofacets_[compType].append(reqCofacets);
        }
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }
}

void Assembler::configureVector(Vector<double>& b) const 
{
  TimeMonitor timer(configTimer());
  Array<VectorSpace<double> > vs(eqn_->numVarBlocks());
  for (unsigned int i=0; i<eqn_->numVarBlocks(); i++)
    {
      vs[i] = rowSpace_[i]->vecSpace();
    }
  VectorSpace<double> rowSpace;

  if ((int) vs.size() > 1)
    {
      rowSpace = TSFExtended::productSpace(vs);
    }
  else
    {
      rowSpace = vs[0];
    }

  b = rowSpace.createMember();

  if (rowSpace.numBlocks() > 1)
    {
      /* configure the blocks */
      Vector<double> vecBlock;
      for (int br=0; br<rowSpace.numBlocks(); br++)
        {
          configureVectorBlock(br, vecBlock);
          b.setBlock(br, vecBlock);
        }
    }
  else
    {
      /* nothing to do here except check that the vector is loadable */
      TSFExtended::LoadableVector<double>* lv 
        = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());
  
      TEST_FOR_EXCEPTION(lv == 0, RuntimeError,
                         "vector is not loadable in Assembler::configureVector()");
    }
  
  vecNeedsConfiguration_ = false;
}

void Assembler::configureVectorBlock(int br, Vector<double>& b) const 
{
  VectorSpace<double> vecSpace = rowSpace_[br]->vecSpace();

  b = vecSpace.createMember();
  
  TSFExtended::LoadableVector<double>* lv 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());
  
  TEST_FOR_EXCEPTION(lv == 0, RuntimeError,
                     "vector block is not loadable "
                     "in Assembler::configureVectorBlock()");
}


void Assembler::configureMatrix(LinearOperator<double>& A,
                                Vector<double>& b) const
{
  TimeMonitor timer(configTimer());
  int nRowBlocks = rowMap_.size();
  int nColBlocks = colMap_.size();
  Array<Array<int> > isNonzero = findNonzeroBlocks();

  if (nRowBlocks==1 && nColBlocks==1)
    {
      configureMatrixBlock(0,0,A);
    }
  else
    {
      A = new BlockOperator<double>(solnVecSpace(), rowVecSpace());
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
    }
  
  configureVector(b);

  matNeedsConfiguration_ = false;
}

void Assembler::configureMatrixBlock(int br, int bc,
                                     LinearOperator<double>& A) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());
  
  SUNDANCE_VERB_LOW(tab << "Assembler: num rows = " << rowMap()[br]->numDOFs());
  
  SUNDANCE_VERB_LOW(tab << "Assembler: num cols = " << colMap()[bc]->numDOFs());

  VectorSpace<double> rowSpace = rowSpace_[br]->vecSpace();
  VectorSpace<double> colSpace = colSpace_[bc]->vecSpace();

  RefCountPtr<MatrixFactory<double> > matFactory 
    = rowVecType_[br].createMatrixFactory(colSpace, rowSpace);

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(matFactory.get());

  CollectivelyConfigurableMatrixFactory* ccmf 
    = dynamic_cast<CollectivelyConfigurableMatrixFactory*>(matFactory.get());

  TEST_FOR_EXCEPTION(ccmf==0 && icmf==0, RuntimeError,
                     "Neither incremental nor collective matrix structuring "
                     "appears to be available");


  /* If collective structuring is the user preference, or if incremental
   * structuring is not supported, do collective structuring */
  if ((icmf==0 || !matrixEliminatesRepeatedCols()) && ccmf != 0)
    {
      SUNDANCE_VERB_MEDIUM(tab << "Assembler: doing collective matrix structuring...");
      Array<int> graphData;
      Array<int> nnzPerRow;
      Array<int> rowPtrs;

      getGraph(br, bc, graphData, rowPtrs, nnzPerRow);
      ccmf->configure(lowestRow_[br], rowPtrs, nnzPerRow, graphData);
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tab << "Assembler: doing incremental matrix structuring...");
      incrementalGetGraph(br, bc, icmf);
      {
        TimeMonitor timer1(matFinalizeTimer());
        icmf->finalize();
      }
    }
  
  SUNDANCE_VERB_MEDIUM(tab << "Assembler: done");

  SUNDANCE_VERB_MEDIUM(tab << "Assembler: constructing matrix...");
  {
    TimeMonitor timer1(matAllocTimer());
    A = matFactory->createMatrix();
  }

  SUNDANCE_VERB_MEDIUM(tab << "...done");
}


/* ------------  assemble both the vector and the matrix  ------------- */

void Assembler::assemble(LinearOperator<double>& A,
                         Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  numAssembleCalls()++;

  TEST_FOR_EXCEPTION(!contexts_.containsKey(MatrixAndVector),
                     RuntimeError,
                     "Assembler::assemble(A, b) called for an assembler that "
                     "does not support matrix/vector assembly");

  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  SUNDANCE_VERB_LOW(tab << "Assembling matrix and vector"); 

  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());
  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;

  if (matNeedsConfiguration_)
    {
      configureMatrix(A, b);
    }
  
  int numRowBlocks = A.range().numBlocks();
  int numColBlocks = A.domain().numBlocks();

  RefCountPtr<Array<Array<Array<int> > > > testLocalDOFs 
    = rcp(new Array<Array<Array<int> > >(numRowBlocks));

  RefCountPtr<Array<Array<Array<int> > > > unkLocalDOFs
    = rcp(new Array<Array<Array<int> > >(numColBlocks));

  Array<TSFExtended::LoadableVector<double>* > vec(numRowBlocks);
  Array<Array<TSFExtended::LoadableMatrix<double>* > > mat(numRowBlocks, numColBlocks);

  for (int br=0; br<numRowBlocks; br++)
    {
      Vector<double> vecBlock = b.getBlock(br);
      vec[br] = dynamic_cast<TSFExtended::LoadableVector<double>* >(vecBlock.ptr().get());
      TEST_FOR_EXCEPTION(vec[br]==0, RuntimeError,
                         "vector block " << br 
                         << " is not loadable in Assembler::assemble()");
      vecBlock.zero();
      for (int bc=0; bc<numColBlocks; bc++)
        {
          LinearOperator<double> matBlock = A.getBlock(br, bc);
          if (matBlock.ptr().get() == 0) continue;
          mat[br][bc] 
            = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(matBlock.ptr().get());
          TEST_FOR_EXCEPTION(mat[br][bc]==0, RuntimeError,
                             "matrix block (" << br << ", " << bc 
                             << ") is not loadable in Assembler::assemble()");
          mat[br][bc]->zero();
        }
      if (verbosity() > VerbHigh)
        {
          Tabs tab1;
          cerr << tab1 << "map" << endl;
          rowMap_[br]->print(cerr);
          cerr << tab1 << "BC row flags " << endl;
          for (unsigned int i=0; i<isBCRow_[br]->size(); i++) 
            {
              cerr << tab1 << i << " " << (*(isBCRow_[br]))[i] << endl;
            }
        }
    }

  const Array<EvalContext>& contexts = contexts_.get(MatrixAndVector);
  const Array<Array<IntegralGroup> >& groups = groups_.get(MatrixAndVector);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(MatrixAndVector);

  for (unsigned int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_VERB_MEDIUM(endl 
                           << "================================================="
                           << endl << tab0 << " doing subregion=" 
                           << rqc_[r]);     


      SUNDANCE_VERB_MEDIUM(tab0 << "expr is " << evalExprs[r]->toString());
      SUNDANCE_VERB_MEDIUM(tab0 << "isBC= " << isBCRqc_[r]);

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts_.get(MatrixAndVector)[r]);

      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);
      mediators_[r]->setCellType(cellType);      

      SUNDANCE_VERB_MEDIUM(tab0 << "cell type = " << cellType);

      const Evaluator* evaluator 
        = evalExprs[r]->evaluator(contexts[r]).get();

      /* do the cells in batches of the work set size */

      CellIterator iter=cells.begin();
      int workSetCounter = 0;

      while (iter != cells.end())
        {
          Tabs tab1;
          /* build up the work set */
          workSet->resize(0);
          for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }
          SUNDANCE_VERB_MEDIUM(
                               tab1 << "doing work set " << workSetCounter
                               << " consisting of " << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);

          workSetCounter++;

          bool useMaximalCellsForTransformations 
            = rqcRequiresMaximalCofacets_.get(MatrixAndVector)[r];
          mediators_[r]->setCellBatch(useMaximalCellsForTransformations, workSet);

          const CellJacobianBatch& JVol = mediators_[r]->JVol();
          const CellJacobianBatch& JTrans = mediators_[r]->JTrans();
          const Array<int>& facetIndices = mediators_[r]->facetIndices();

          evaluator->resetNumCalls();
          evalExprs[r]->evaluate(*evalMgr_, constantCoeffs, vectorCoeffs);

          if (verbosity() > VerbHigh)
            {
              Tabs tab2;
              cerr << tab2 << "evaluation results: " << endl;
              const EvalContext& context = contexts[r];
              const RefCountPtr<SparsitySuperset>& sparsity 
                = evalExprs[r]->sparsitySuperset(context);
              sparsity->print(cerr, vectorCoeffs, constantCoeffs);
            }

          Array<Array<int> > nTestNodes(numRowBlocks);
          Array<Array<int> > nUnkNodes(numColBlocks);
          
          for (int br=0; br<numRowBlocks; br++)
            {   
              rowMap_[br]->getDOFsForCellBatch(cellDim, *workSet, 
                                               (*testLocalDOFs)[br], nTestNodes[br]);
            }
          SUNDANCE_VERB_EXTREME(tab1 << "local DOF values " << *testLocalDOFs);
          for (int bc=0; bc<numColBlocks; bc++)
            {
              if (bc < numRowBlocks && rowMap_[bc].get()==colMap_[bc].get())
                {
                  (*unkLocalDOFs)[bc] = (*testLocalDOFs)[bc];
                  nUnkNodes[bc] = nTestNodes[bc];
                }
              else
                {
                  colMap_[bc]->getDOFsForCellBatch(cellDim, *workSet,
                                                   (*unkLocalDOFs)[bc],
                                                   nUnkNodes[bc]);
                }
            }
          
          ElementIntegral::invalidateTransformationMatrices();

          for (unsigned int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(JTrans, JVol, facetIndices, vectorCoeffs,
                                  constantCoeffs, 
                                  localValues)) continue;

              if (verbosity() > VerbHigh)
                {
                  cerr << endl << endl 
                       << "--------------- doing integral group " << g << endl;
                  cerr << "num test DOFs = " << testLocalDOFs->size() << endl;
                  cerr << "num unk DOFs = " << unkLocalDOFs->size() << endl;
                  cerr << "num entries = " << localValues->size() << endl;
                  cerr << "values = " << *localValues << endl;
                }
              if (group.isTwoForm())
                {
                  insertLocalMatrixBatch(cellDim, *workSet, isBCRqc_[r],
                                         *testLocalDOFs,
                                         *unkLocalDOFs,
                                         nTestNodes,
                                         nUnkNodes,
                                         group.testID(), group.testBlock(),
                                         group.unkID(), group.unkBlock(),
                                         *localValues, mat);
                }
              else
                {
                  insertLocalVectorBatch(cellDim, *workSet, isBCRqc_[r], 
                                         *testLocalDOFs,
                                         nTestNodes,
                                         group.testID(),group.testBlock(),
                                         *localValues, vec);
                }
            }
        }
    }

  SUNDANCE_VERB_LOW(tab << "Assembler: done assembling matrix & vector");

  if (verbosity() > VerbHigh)
    {
      cerr << "matrix = " << endl;
      A.print(cerr);
      cerr << "vector = " << endl;
      b.print(cerr);
    }

}



/* ------------  assemble the vector alone  ------------- */

void Assembler::assemble(Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  numAssembleCalls()++;
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  TEST_FOR_EXCEPTION(!contexts_.containsKey(VectorOnly),
                     RuntimeError,
                     "Assembler::assemble(b) called for an assembler that "
                     "does not support vector-only assembly");

  SUNDANCE_VERB_LOW(tab << "Assembling vector"); 

  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;

  if (vecNeedsConfiguration_)
    {
      configureVector(b);
    }

  
  int numRowBlocks = b.space().numBlocks();

  RefCountPtr<Array<Array<Array<int> > > > testLocalDOFs 
    = rcp(new Array<Array<Array<int> > >(numRowBlocks));

  Array<TSFExtended::LoadableVector<double>* > vec(numRowBlocks);

  for (int br=0; br<numRowBlocks; br++)
    {
      Vector<double> vecBlock = b.getBlock(br);
      vec[br] = dynamic_cast<TSFExtended::LoadableVector<double>* >(vecBlock.ptr().get());
      TEST_FOR_EXCEPTION(vec[br]==0, RuntimeError,
                         "vector block " << br 
                         << " is not loadable in Assembler::assemble()");
      vecBlock.zero();
    }


  const Array<EvalContext>& contexts = contexts_.get(VectorOnly);
  const Array<Array<IntegralGroup> >& groups = groups_.get(VectorOnly);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(VectorOnly);

  for (unsigned int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;

      SUNDANCE_VERB_MEDIUM(tab0 << "doing subregion=" 
                           << rqc_[r]);     


      
      SUNDANCE_VERB_MEDIUM(tab0 << "expr is " << evalExprs[r]->toString());

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts[r]);

      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);
      mediators_[r]->setCellType(cellType);      


      SUNDANCE_VERB_MEDIUM(tab0 << "cell type = " << cellType);

      const Evaluator* evaluator 
        = evalExprs[r]->evaluator(contexts[r]).get();
      /* do the cells in batches of the work set size */

      CellIterator iter=cells.begin();
      int workSetCounter = 0;

      while (iter != cells.end())
        {
          Tabs tab1;
          /* build up the work set */
          workSet->resize(0);
          for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);


          workSetCounter++;
          bool useMaximalCellsForTransformations 
            = rqcRequiresMaximalCofacets_.get(VectorOnly)[r];
          mediators_[r]->setCellBatch(useMaximalCellsForTransformations, workSet);

          const CellJacobianBatch& JVol = mediators_[r]->JVol();
          const CellJacobianBatch& JTrans = mediators_[r]->JTrans();
          const Array<int>& facetIndices = mediators_[r]->facetIndices();

          evaluator->resetNumCalls();
          //          mesh_.getJacobians(cellDim, *workSet, *J);
          evalExprs[r]->evaluate(*evalMgr_, constantCoeffs, vectorCoeffs);

          if (verbosity() > VerbHigh)
            {
              Tabs tab2;
              cerr << tab2 << " ----------- evaluation results: ------" << endl;
              cerr << tab2 << "expr=" << evalExprs[r]->toString() << endl;
              const EvalContext& context = contexts[r];
              const RefCountPtr<SparsitySuperset>& sparsity 
                = evalExprs[r]->sparsitySuperset(context);
              sparsity->print(cerr, vectorCoeffs, constantCoeffs);
            }

          Array<Array<int> > nTestNodes(numRowBlocks);
          for (int br=0; br<numRowBlocks; br++)
            {
              rowMap_[br]->getDOFsForCellBatch(cellDim, *workSet, (*testLocalDOFs)[br],
                                               nTestNodes[br]);
            }

          ElementIntegral::invalidateTransformationMatrices();
          for (unsigned int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(JTrans, JVol, facetIndices, vectorCoeffs, 
                                  constantCoeffs, 
                                  localValues)) 
                {
                  continue;
                }

              if (verbosity() > VerbHigh)
                {
                  cerr << endl << endl 
                       << "--------------- doing integral group " << g << endl;
                  cerr << "num test DOFs = " << testLocalDOFs->size() << endl;
                  cerr << "num entries = " << localValues->size() << endl;
                  cerr << "values = " << *localValues << endl;
                }

              insertLocalVectorBatch(cellDim, *workSet, isBCRqc_[r], 
                                     *testLocalDOFs,
                                     nTestNodes,
                                     group.testID(), group.testBlock(), *localValues, vec);
            }
        }
    }

  SUNDANCE_VERB_LOW(tab << "Assembler: done assembling vector");
  if (verbosity() > VerbHigh)
    {
      cerr << "vector = " << endl;
      b.print(cerr);
    }
}


/* ------------  evaluate a functional and its gradient ---- */

void Assembler::evaluate(double& value, Vector<double>& gradient) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  numAssembleCalls()++;
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  TEST_FOR_EXCEPTION(!contexts_.containsKey(FunctionalAndGradient),
                     RuntimeError,
                     "Assembler::evaluate(f,df) called for an assembler that "
                     "does not support value/gradient assembly");

  SUNDANCE_VERB_LOW("-----------------------------------------------------------------"); 
  SUNDANCE_VERB_LOW("------- Computing functional and gradient "); 
  SUNDANCE_VERB_LOW("-----------------------------------------------------------------"); 

  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;

  if (vecNeedsConfiguration_)
    {
      configureVector(gradient);
    }

  
  int numRowBlocks = gradient.space().numBlocks();

  RefCountPtr<Array<Array<Array<int> > > > testLocalDOFs 
    = rcp(new Array<Array<Array<int> > >(numRowBlocks));

  Array<TSFExtended::LoadableVector<double>* > vec(numRowBlocks);

  for (int br=0; br<numRowBlocks; br++)
    {
      Vector<double> vecBlock = gradient.getBlock(br);
      vec[br] = dynamic_cast<TSFExtended::LoadableVector<double>* >(vecBlock.ptr().get());
      TEST_FOR_EXCEPTION(vec[br]==0, RuntimeError,
                         "vector block " << br 
                         << " is not loadable in Assembler::assemble()");
      vecBlock.zero();
    }

  double localSum = 0.0;

  const Array<EvalContext>& contexts = contexts_.get(FunctionalAndGradient);
  const Array<Array<IntegralGroup> >& groups 
    = groups_.get(FunctionalAndGradient);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(FunctionalAndGradient);

  for (unsigned int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;

      SUNDANCE_VERB_MEDIUM(tab0 << "doing subregion=" 
                           << rqc_[r]);     


      
      SUNDANCE_VERB_MEDIUM(tab0 << "expr is " << evalExprs[r]->toString());

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts[r]);

      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);
      mediators_[r]->setCellType(cellType);      


      SUNDANCE_VERB_MEDIUM(tab0 << "cell type = " << cellType);

      const Evaluator* evaluator 
        = evalExprs[r]->evaluator(contexts[r]).get();
      /* do the cells in batches of the work set size */

      CellIterator iter=cells.begin();
      int workSetCounter = 0;

      while (iter != cells.end())
        {
          Tabs tab1;
          /* build up the work set */
          workSet->resize(0);
          for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);


          workSetCounter++;

          bool useMaximalCellsForTransformations 
            = rqcRequiresMaximalCofacets_.get(FunctionalAndGradient)[r];
          mediators_[r]->setCellBatch(useMaximalCellsForTransformations, workSet);

          const CellJacobianBatch& JVol = mediators_[r]->JVol();
          const CellJacobianBatch& JTrans = mediators_[r]->JTrans();
          const Array<int>& facetIndices = mediators_[r]->facetIndices();


          evaluator->resetNumCalls();
          //          mesh_.getJacobians(cellDim, *workSet, *J);
          evalExprs[r]->evaluate(*evalMgr_, constantCoeffs, vectorCoeffs);

          if (verbosity() > VerbHigh)
            {
              Tabs tab2;
              cerr << tab2 << " ----------- evaluation results: ------" << endl;
              cerr << tab2 << "expr=" << evalExprs[r]->toString() << endl;
              const EvalContext& context = contexts[r];
              const RefCountPtr<SparsitySuperset>& sparsity 
                = evalExprs[r]->sparsitySuperset(context);
              sparsity->print(cerr, vectorCoeffs, constantCoeffs);
            }

          Array<Array<int> > nTestNodes(numRowBlocks);
          for (int br=0; br<numRowBlocks; br++)
            {
              rowMap_[br]->getDOFsForCellBatch(cellDim, *workSet, (*testLocalDOFs)[br],
                                               nTestNodes[br]);
            }

          ElementIntegral::invalidateTransformationMatrices();

          for (unsigned int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(JTrans, JVol, facetIndices, vectorCoeffs, 
                                  constantCoeffs, 
                                  localValues)) 
                {
                  continue;
                }

              if (verbosity() > VerbHigh)
                {
                  cerr << endl << endl 
                       << "--------------- doing integral group " << g << endl;
                  cerr << "num test DOFs = " << testLocalDOFs->size() << endl;
                  cerr << "num entries = " << localValues->size() << endl;
                  cerr << "values = " << *localValues << endl;
                }
              
              if (group.isOneForm())
                {
                  insertLocalVectorBatch(cellDim, *workSet, isBCRqc_[r], 
                                         *testLocalDOFs,
                                         nTestNodes,
                                         group.testID(), group.testBlock(), 
                                         *localValues, vec);
                }
              else
                {
                  localSum += (*localValues)[0];
                }
            }
        }
    }

  value = localSum;

  mesh_.comm().allReduce((void*) &localSum, (void*) &value, 1, 
                         MPIComm::DOUBLE, MPIComm::SUM);

  SUNDANCE_VERB_LOW(tab << "Assembler: done computing functional and its gradient");
  if (verbosity() > VerbHigh)
    {
      cerr << "vector = " << endl;
      gradient.print(cerr);
    }
}




/* ------------  evaluate a functional ---- */

void Assembler::evaluate(double& value) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());
  numAssembleCalls()++;
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  TEST_FOR_EXCEPTION(!contexts_.containsKey(FunctionalOnly),
                     RuntimeError,
                     "Assembler::evaluate(f) called for an assembler that "
                     "does not support functional evaluation");


  SUNDANCE_VERB_LOW("-----------------------------------------------------------------"); 
  SUNDANCE_VERB_LOW("------- Computing functional"); 
  SUNDANCE_VERB_LOW("-----------------------------------------------------------------"); 


  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;

  double localSum = 0.0;

  const Array<EvalContext>& contexts = contexts_.get(FunctionalOnly);
  const Array<Array<IntegralGroup> >& groups 
    = groups_.get(FunctionalOnly);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(FunctionalOnly);

  for (unsigned int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;

      SUNDANCE_VERB_MEDIUM(tab0 << "doing subregion=" 
                           << rqc_[r]);     


      
      SUNDANCE_VERB_MEDIUM(tab0 << "expr is " << evalExprs[r]->toString());

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts[r]);

      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);
      mediators_[r]->setCellType(cellType);      


      SUNDANCE_VERB_MEDIUM(tab0 << "cell type = " << cellType);

      const Evaluator* evaluator 
        = evalExprs[r]->evaluator(contexts[r]).get();
      /* do the cells in batches of the work set size */

      CellIterator iter=cells.begin();
      int workSetCounter = 0;

      while (iter != cells.end())
        {
          Tabs tab1;
          /* build up the work set */
          workSet->resize(0);
          for (unsigned int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);



          bool useMaximalCellsForTransformations 
            = rqcRequiresMaximalCofacets_.get(FunctionalOnly)[r];
          mediators_[r]->setCellBatch(useMaximalCellsForTransformations, workSet);

          const CellJacobianBatch& JVol = mediators_[r]->JVol();
          const CellJacobianBatch& JTrans = mediators_[r]->JTrans();
          const Array<int>& facetIndices = mediators_[r]->facetIndices();


          evaluator->resetNumCalls();
          evalExprs[r]->evaluate(*evalMgr_, constantCoeffs, vectorCoeffs);

          if (verbosity() > VerbHigh)
            {
              Tabs tab2;
              cerr << tab2 << " ----------- evaluation results: ------" << endl;
              cerr << tab2 << "expr=" << evalExprs[r]->toString() << endl;
              const EvalContext& context = contexts[r];
              const RefCountPtr<SparsitySuperset>& sparsity 
                = evalExprs[r]->sparsitySuperset(context);
              sparsity->print(cerr, vectorCoeffs, constantCoeffs);
            }

          ElementIntegral::invalidateTransformationMatrices();
          for (unsigned int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(JTrans, JVol, facetIndices, vectorCoeffs, 
                                  constantCoeffs, 
                                  localValues)) 
                {
                  continue;
                }
              SUNDANCE_VERB_HIGH(tab1 << "contribution from work set "
                                 << workSetCounter << " is " 
                                 << (*localValues)[0]);
              localSum += (*localValues)[0];
            }
          workSetCounter++;
        }
    }

  value = localSum;

  mesh_.comm().allReduce((void*) &localSum, (void*) &value, 1, 
                         MPIComm::DOUBLE, MPIComm::SUM);
  SUNDANCE_VERB_LOW(tab << "Assembler: done computing functional");
}








/* ------------  insert elements into the matrix  ------------- */


void Assembler::insertLocalMatrixBatch(int cellDim, 
                                       const Array<int>& workSet, 
                                       bool isBCRqc,
                                       const Array<Array<Array<int> > >& testIndices,
                                       const Array<Array<Array<int> > >& unkIndices,
                                       const Array<Array<int> >& nTestNodes, 
                                       const Array<Array<int> >& nUnkNodes,
                                       const Array<int>& testID, 
                                       const Array<int>& testBlock, 
                                       const Array<int>& unkID,
                                       const Array<int>& unkBlock,
                                       const Array<double>& localValues, 
                                       Array<Array<LoadableMatrix<double>* > >& mat) const 
{
  Tabs tab;
  TimeMonitor timer(matInsertTimer());

  SUNDANCE_VERB_HIGH(tab << "inserting local matrix values...");

  static Array<int> skipRow;
  static Array<int> rows;
  static Array<int> cols;

  int nCells = workSet.size();


  if (verbosity() > VerbHigh)
    {
      cerr << "isBC " << isBCRqc << endl;
      cerr << "num test nodes = " << nTestNodes << endl;
      cerr << "num unk nodes = " << nUnkNodes << endl;
      cerr << "num cells = " << nCells << endl;
      cerr << "testID = " << testID << endl;
      cerr << "unkID = " << unkID << endl;
    }
  for (unsigned int t=0; t<testID.size(); t++)
    {
      Tabs tab1;
      SUNDANCE_VERB_EXTREME(tab1 << "test ID = " << testID[t]);
      SUNDANCE_VERB_EXTREME(tab1 << "is BC eqn = " << isBCRqc);
      int br = testBlock[t];
      const RefCountPtr<DOFMapBase>& rowMap = rowMap_[br];
      int highestIndex = lowestRow_[br] + rowMap->numLocalDOFs();
      int lowestLocalRow = rowMap->lowestLocalDOF();
      int testChunk = rowMap->chunkForFuncID(testID[t]);
      int testFuncIndex = rowMap->indexForFuncID(testID[t]);
      const Array<int>& testDOFs = testIndices[br][testChunk];
      SUNDANCE_VERB_EXTREME(tab1 << "test DOFs = " << testDOFs);
      int nTestFuncs = rowMap->nFuncs(testChunk);
      int numTestNodes = nTestNodes[br][testChunk];
      int numRows = nCells * numTestNodes;
      const Array<int>& isBCRow = *(isBCRow_[br]);
      SUNDANCE_VERB_EXTREME(tab1 << "isBCRow = " << isBCRow);
      rows.resize(numRows);
      skipRow.resize(numRows);
      int r=0;
      for (int c=0; c<nCells; c++)
        {
          for (int n=0; n<numTestNodes; n++, r++)
            {
              int row = testDOFs[(c*nTestFuncs + testFuncIndex)*numTestNodes + n];
              rows[r] = row;
              int localRow = rows[r]-lowestLocalRow;
              skipRow[r] = row < lowestLocalRow || row >= highestIndex
                || (isBCRqc && !isBCRow[localRow])
                || (!isBCRqc && isBCRow[localRow]);
            }
        }

      if (verbosity() > VerbHigh)
        {
          Tabs tab2;
          Set<int> activeRows;
          Set<int> skippedRows;
          for (int i=0; i<numRows; i++) 
            {
              if (skipRow[i]) skippedRows.put(rows[i]);
              else activeRows.put(rows[i]);
            }
          cerr << tab2 << "active rows = " << activeRows << endl;
          cerr << tab2 << "skipped rows = " << skippedRows << endl;
        }

      for (unsigned int u=0; u<unkID.size(); u++)
        {      
          int bc = unkBlock[u];
          const RefCountPtr<DOFMapBase>& colMap = colMap_[bc];
          int unkChunk = colMap->chunkForFuncID(unkID[u]);
          int unkFuncIndex = colMap->indexForFuncID(unkID[u]);
          const Array<int>& unkDOFs = unkIndices[bc][unkChunk];
          int nUnkFuncs = colMap->nFuncs(unkChunk);
          int numUnkNodes = nUnkNodes[bc][unkChunk];
          cols.resize(nCells*numUnkNodes);
          int j=0;
          for (int c=0; c<nCells; c++)
            {
              for (int n=0; n<numUnkNodes; n++, j++)
                {
                  cols[j] = unkDOFs[(c*nUnkFuncs + unkFuncIndex)*numUnkNodes + n];
                }
            }
          
          mat[br][bc]->addToElementBatch(numRows,
                                         numTestNodes,
                                         &(rows[0]),
                                         numUnkNodes,
                                         &(cols[0]),
                                         &(localValues[0]),
                                         &(skipRow[0]));
        }
    }
}







/* ------------  insert elements into the vector  ------------- */

void Assembler
::insertLocalVectorBatch(int cellDim, 
                         const Array<int>& workSet, 
                         bool isBCRqc,
                         const Array<Array<Array<int> > >& testIndices,
                         const Array<Array<int> >& nTestNodes, 
                         const Array<int>& testID,  
                         const Array<int>& testBlock, 
                         const Array<double>& localValues, 
                         Array<TSFExtended::LoadableVector<double>* >& vec) const 
{
  TimeMonitor timer(vecInsertTimer());
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "inserting local vector values...");
  SUNDANCE_VERB_EXTREME(tab << "values are " << localValues);


  int nCells = workSet.size();

  for (unsigned int i=0; i<testID.size(); i++)
    {
      Tabs tab1;
      SUNDANCE_VERB_EXTREME(tab1 << "test ID = " << testID[i]);
      SUNDANCE_VERB_EXTREME(tab1 << "is BC eqn = " << isBCRqc);
      int br = testBlock[i];
      const RefCountPtr<DOFMapBase>& rowMap = rowMap_[br];
      int lowestLocalRow = rowMap->lowestLocalDOF();
      int chunk = rowMap->chunkForFuncID(testID[i]);
      int funcIndex = rowMap->indexForFuncID(testID[i]);
      const Array<int>& dofs = testIndices[br][chunk];
      int nFuncs = rowMap->nFuncs(chunk);
      int nNodes = nTestNodes[br][chunk];
      const Array<int>& isBCRow = *(isBCRow_[br]);
      int r=0;
      TSFExtended::LoadableVector<double>* vecBlock = vec[br];

      for (int c=0; c<nCells; c++)
        {
          for (int n=0; n<nNodes; n++, r++)
            {
              int rowIndex = dofs[(c*nFuncs + funcIndex)*nNodes + n];
              int localRowIndex = rowIndex - lowestLocalRow;
              if (!(rowMap->isLocalDOF(rowIndex))
                  || isBCRqc!=isBCRow[localRowIndex]) continue;
              {
                vecBlock->addToElement(rowIndex, localValues[r]);
              }
            }
        }
    }
  SUNDANCE_VERB_HIGH(tab << "...done");
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
        SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                     tab0 << "cell set " << domain
                     << " isBCRegion=" << eqn_->isBCRegion(d));
        unsigned int dim = domain.dimension(mesh_);
        CellSet cells = domain.getCells(mesh_);

        RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
        if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

        SUNDANCE_OUT(this->verbosity() > VerbMedium && pairs.get() != 0, 
                     tab0 << "non-BC pairs = "
                     << *pairs);
       
        RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
        if (eqn_->isBCRegion(d))
          {
            if (eqn_->hasBCVarUnkPairs(domain)) 
              {
                bcPairs = eqn_->bcVarUnkPairs(domain);
                SUNDANCE_OUT(this->verbosity() > VerbMedium, tab0 << "BC pairs = "
                             << *bcPairs);
              }
          }
        Array<Set<int> > unksForTestsSet(eqn_->numVars(bc));
        Array<Set<int> > bcUnksForTestsSet(eqn_->numVars(bc));

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

        int nt = eqn_->numVars(br);
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

            rowMap_[br]->getDOFsForCellBatch(dim, *workSet, *testLocalDOFs,
                                             numTestNodes);
            if (rowMap_[br].get()==colMap_[bc].get())
              {
                unkLocalDOFs = testLocalDOFs;
                numUnkNodes = numTestNodes;
              }
            else
              {
                colMap_[br]->getDOFsForCellBatch(dim, *workSet, 
                                                 *unkLocalDOFs, numUnkNodes);
              }

            if (pairs.get() != 0)
              {
                for (int c=0; c<nCells; c++)
                  {
                    for (int t=0; t<nt; t++)
                      {
                        int tChunk = rowMap_[br]->chunkForFuncID(t);
                        int nTestFuncs = rowMap_[br]->nFuncs(tChunk);
                        int testFuncIndex = rowMap_[br]->indexForFuncID(t);
                        int nTestNodes = numTestNodes[tChunk];
                        const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
                        for (unsigned int uit=0; uit<unksForTests[t].size(); uit++)
                          {
                            Tabs tab2;
                            int u = unksForTests[t][uit];
                            int uChunk = colMap_[bc]->chunkForFuncID(u);
                            int nUnkFuncs = colMap_[bc]->nFuncs(uChunk);
                            int unkFuncIndex = colMap_[bc]->indexForFuncID(u);
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
                        int tChunk = rowMap_[br]->chunkForFuncID(t);
                        int nTestFuncs = rowMap_[br]->nFuncs(tChunk);
                        int testFuncIndex = rowMap_[br]->indexForFuncID(t);
                        int nTestNodes = numTestNodes[tChunk];
                        const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
                        for (unsigned int uit=0; uit<bcUnksForTests[t].size(); uit++)
                          {
                            Tabs tab2;
                            int u = bcUnksForTests[t][uit];
                            int uChunk = colMap_[bc]->chunkForFuncID(u);
                            int nUnkFuncs = colMap_[bc]->nFuncs(uChunk);
                            int unkFuncIndex = colMap_[bc]->indexForFuncID(u);
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
      SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                   tab0 << "cell set " << domain
                   << " isBCRegion=" << eqn_->isBCRegion(d));
      unsigned int dim = domain.dimension(mesh_);
      CellSet cells = domain.getCells(mesh_);

      RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
      if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

      SUNDANCE_OUT(this->verbosity() > VerbMedium && pairs.get() != 0, 
                   tab0 << "non-BC pairs = "
                   << *pairs);
       
      RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
      if (eqn_->isBCRegion(d))
        {
          if (eqn_->hasBCVarUnkPairs(domain)) 
            {
              bcPairs = eqn_->bcVarUnkPairs(domain);
              SUNDANCE_OUT(this->verbosity() > VerbMedium, tab0 << "BC pairs = "
                           << *bcPairs);
            }
        }
      Array<Set<int> > unksForTestsSet(eqn_->numVars(br));
      Array<Set<int> > bcUnksForTestsSet(eqn_->numVars(bc));

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

      int nt = eqn_->numVars(br);
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

          rowMap_[br]->getDOFsForCellBatch(dim, *workSet, *testLocalDOFs,
                                       numTestNodes);
          if (rowMap_[br].get()==colMap_[bc].get())
            {
              unkLocalDOFs = testLocalDOFs;
              numUnkNodes = numTestNodes;
            }
          else
            {
              colMap_[bc]->getDOFsForCellBatch(dim, *workSet, 
                                           *unkLocalDOFs, numUnkNodes);
            }

          
          if (pairs.get() != 0)
            {
              for (int c=0; c<nCells; c++)
                {
                  for (int t=0; t<nt; t++)
                    {
                      int tChunk = rowMap_[br]->chunkForFuncID(t);
                      int nTestFuncs = rowMap_[br]->nFuncs(tChunk);
                      int testFuncIndex = rowMap_[br]->indexForFuncID(t);
                      int nTestNodes = numTestNodes[tChunk];
                      const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
                      for (unsigned int uit=0; uit<unksForTests[t].size(); uit++)
                        {
                          Tabs tab2;
                          int u = unksForTests[t][uit];
                          int uChunk = colMap_[bc]->chunkForFuncID(u);
                          int nUnkFuncs = colMap_[bc]->nFuncs(uChunk);
                          int unkFuncIndex = colMap_[bc]->indexForFuncID(u);
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
                  for (int t=0; t<nt; t++)
                    {
                      int tChunk = rowMap_[br]->chunkForFuncID(t);
                      int nTestFuncs = rowMap_[br]->nFuncs(tChunk);
                      int testFuncIndex = rowMap_[br]->indexForFuncID(t);
                      int nTestNodes = numTestNodes[tChunk];
                      const Array<int>& testDOFs = (*testLocalDOFs)[tChunk];
                      for (unsigned int uit=0; uit<bcUnksForTests[t].size(); uit++)
                        {
                          Tabs tab2;
                          int u = bcUnksForTests[t][uit];
                          int uChunk = colMap_[bc]->chunkForFuncID(u);
                          int nUnkFuncs = colMap_[bc]->nFuncs(uChunk);
                          int unkFuncIndex = colMap_[bc]->indexForFuncID(u);
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
  Array<Array<int> > rtn(eqn_->numVarBlocks(), eqn_->numUnkBlocks());
  for (unsigned int br=0; br<rtn.size(); br++)
    {
      for (unsigned int bc=0; bc<rtn[br].size(); bc++)
        {
          rtn[br][bc] = 0 ;
        }
    }

  for (unsigned int d=0; d<eqn_->numRegions(); d++)
    {
      Tabs tab0;
      CellFilter domain = eqn_->region(d);
      SUNDANCE_OUT(this->verbosity() > VerbMedium, 
                   tab0 << "cell set " << domain
                   << " isBCRegion=" << eqn_->isBCRegion(d));

      RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
      if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

      SUNDANCE_OUT(this->verbosity() > VerbMedium && pairs.get() != 0, 
                   tab0 << "non-BC pairs = "
                   << *pairs);
       
      RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
      if (eqn_->isBCRegion(d))
        {
          if (eqn_->hasBCVarUnkPairs(domain)) 
            {
              bcPairs = eqn_->bcVarUnkPairs(domain);
              SUNDANCE_OUT(this->verbosity() > VerbMedium, tab0 << "BC pairs = "
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



