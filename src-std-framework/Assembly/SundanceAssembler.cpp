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

static Time& graphFlatteningTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("tmp graph flattening"); 
  return *rtn;
}


Assembler
::Assembler(const Mesh& mesh, 
            const RefCountPtr<EquationSet>& eqn,
            const VectorType<double>& vectorType,
            const VerbositySetting& verb)
  : matNeedsConfiguration_(true),
    vecNeedsConfiguration_(true),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    rowSpace_(),
    colSpace_(),
    bcRows_(),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(),
    lowestRow_(),
    vecType_(vectorType)
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
    vecNeedsConfiguration_(true),
    mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    rowSpace_(),
    colSpace_(),
    bcRows_(),
    rqc_(),
    contexts_(),
    isBCRqc_(),
    groups_(),
    mediators_(),
    evalExprs_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(),
    lowestRow_(),
    vecType_()
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
      rowSpace_ = rcp(new DiscreteSpace(mesh, mapBuilder.testBasisArray(), 
                                        rowMap_, vecType_));
      isBCRow_ = mapBuilder.isBCRow();
      lowestRow_ = mapBuilder.rowMap()->lowestLocalDOF();
    }

  if (!eqn->isFunctionalCalculator())
    {
      colMap_ = mapBuilder.colMap();
      colSpace_ = rcp(new DiscreteSpace(mesh, mapBuilder.unkBasisArray(), 
                                        colMap_, vecType_));
      groups_.put(MatrixAndVector, Array<Array<IntegralGroup> >());
      contexts_.put(MatrixAndVector, Array<EvalContext>());
      evalExprs_.put(MatrixAndVector, Array<const EvaluatableExpr*>());
      groups_.put(VectorOnly, Array<Array<IntegralGroup> >());
      contexts_.put(VectorOnly, Array<EvalContext>());
      evalExprs_.put(VectorOnly, Array<const EvaluatableExpr*>());
    }
  else
    {
      groups_.put(FunctionalAndGradient, Array<Array<IntegralGroup> >());
      contexts_.put(FunctionalAndGradient, Array<EvalContext>());
      evalExprs_.put(FunctionalAndGradient, Array<const EvaluatableExpr*>());
      groups_.put(FunctionalOnly, Array<Array<IntegralGroup> >());
      contexts_.put(FunctionalOnly, Array<EvalContext>());
      evalExprs_.put(FunctionalOnly, Array<const EvaluatableExpr*>());
    }

  for (int r=0; r<eqn->regionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];
                         
      rqc_.append(rqc);
      isBCRqc_.append(false);
      const Expr& expr = eqn->expr(rqc);

      SUNDANCE_VERB_HIGH("creating integral groups for rqc=" << rqc << endl
                         << "expr = " << expr);

      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());

      for (Set<ComputationType>::const_iterator 
             i=eqn->computationTypes().begin(); 
           i!=eqn->computationTypes().end();
           i++)
        {
          const ComputationType& compType = *i;
          const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(compType, rqc);
          EvalContext context = eqn->rqcToContext(compType, rqc);
          contexts_[compType].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[compType].append(ee);
          const RefCountPtr<SparsitySuperset>& sparsity 
            = ee->sparsitySuperset(context);
          SUNDANCE_VERB_EXTREME("sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
          groups_[compType].append(groups);
        }
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }

  
  for (int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];
      rqc_.append(rqc);
      isBCRqc_.append(true);
      const Expr& expr = eqn->bcExpr(rqc);

      SUNDANCE_VERB_HIGH("creating integral groups for rqc=" << rqc << endl
                         << "expr = " << expr.toXML().toString());
      
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());

      for (Set<ComputationType>::const_iterator 
             i=eqn->computationTypes().begin(); 
           i!=eqn->computationTypes().end();
           i++)
        {
          const ComputationType& compType = *i;
          const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(compType, rqc);
          EvalContext context = eqn->bcRqcToContext(compType, rqc);
          contexts_[compType].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[compType].append(ee);
          const RefCountPtr<SparsitySuperset>& sparsity 
            = ee->sparsitySuperset(context);
          SUNDANCE_VERB_EXTREME("sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
          groups_[compType].append(groups);
        }
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }
}

void Assembler::configureVector(Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());
  VectorSpace<double> rowSpace = rowSpace_->vecSpace();
  
  b = rowSpace.createMember();

  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::configureVector()");

  vecNeedsConfiguration_ = false;
}



void Assembler::configureMatrix(LinearOperator<double>& A,
                                Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());
  
  SUNDANCE_VERB_LOW(tab << "Assembler: num rows = " << rowMap()->numDOFs());
  
  SUNDANCE_VERB_LOW(tab << "Assembler: num cols = " << colMap()->numDOFs());

  
  SUNDANCE_VERB_MEDIUM(tab << "Assembler: creating row and col spaces...");

  VectorSpace<double> rowSpace = rowSpace_->vecSpace();
  VectorSpace<double> colSpace = colSpace_->vecSpace();

  
  SUNDANCE_VERB_MEDIUM(tab << "...done");




  SUNDANCE_VERB_MEDIUM(tab << "Assembler: creating graph...");

  Array<int> graphData;
  Array<int> nnzPerRow;
  Array<int> rowPtrs;
  getGraph(graphData, rowPtrs, nnzPerRow);

  SUNDANCE_VERB_MEDIUM(tab << "...done");

  SUNDANCE_VERB_MEDIUM(tab << "Assembler: constructing matrix object...");

  A = vecType_.createMatrix(colSpace, rowSpace, &(nnzPerRow[0]));
  //  A = vecType_.createMatrix(colSpace, rowSpace);

  SUNDANCE_VERB_MEDIUM(tab << "...done");
             
  TSFExtended::LoadableMatrix<double>* mat
    = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(A.ptr().get());

  TEST_FOR_EXCEPTION(mat==0, RuntimeError,
                     "matrix is not loadable in Assembler::assemble()");


  SUNDANCE_VERB_MEDIUM(tab << "Assembler: initializing matrix and vector...");


  if (vecNeedsConfiguration_)
    {
      configureVector(b);
    }

  mat->configure(lowestRow_, rowPtrs, nnzPerRow, graphData);

  b.zero();

  SUNDANCE_VERB_MEDIUM(tab << "...done");

  

  matNeedsConfiguration_ = false;
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
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  RefCountPtr<Array<int> > testLocalDOFs 
    = rcp(new Array<int>());

  RefCountPtr<Array<int> > unkLocalDOFs
    = rcp(new Array<int>());

  if (matNeedsConfiguration_)
    {
      configureMatrix(A, b);
    }
  

  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::assemble()");

  TSFExtended::LoadableMatrix<double>* mat
    = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(A.ptr().get());

  TEST_FOR_EXCEPTION(mat==0, RuntimeError,
                     "matrix is not loadable in Assembler::assemble()");

  /* zero out the matrix and vector */
  b.zero();
  mat->zero();

  /* fill loop */

  if (verbosity() > VerbHigh)
    {
      Tabs tab1;
      cerr << tab1 << "map" << endl;
      rowMap_->print(cerr);
      cerr << tab1 << "BC row flags " << endl;
      for (int i=0; i<isBCRow_->size(); i++) 
        {
          cerr << tab1 << i << " " << (*isBCRow_)[i] << endl;
        }
    }

  const Array<EvalContext>& contexts = contexts_.get(MatrixAndVector);
  const Array<Array<IntegralGroup> >& groups = groups_.get(MatrixAndVector);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(MatrixAndVector);

  for (int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_VERB_MEDIUM(tab0 << "doing subregion=" 
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
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }
          SUNDANCE_VERB_MEDIUM(
                               tab1 << "doing work set " << workSetCounter
                               << " consisting of " << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);

          workSetCounter++;

          mediators_[r]->setCellBatch(workSet, J);

          evaluator->resetNumCalls();
          //          mesh_.getJacobians(cellDim, *workSet, *J);
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

          int nTestNodes;
          int nUnkNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, *testLocalDOFs,
                                       nTestNodes);
          SUNDANCE_VERB_EXTREME(tab1 << "local DOF values " << *testLocalDOFs);
          if (rowMap_.get()==colMap_.get())
            {
              unkLocalDOFs = testLocalDOFs;
              nUnkNodes = nTestNodes;
            }
          else
            {
              colMap_->getDOFsForCellBatch(cellDim, *workSet,
                                           *unkLocalDOFs,
                                           nUnkNodes);
            }
          
          ElementIntegral::invalidateTransformationMatrices();

          for (int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(*J, vectorCoeffs,
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
                                         group.nTestNodes(),
                                         group.nUnkNodes(),
                                         group.testID(), group.unkID(), 
                                         *localValues, mat);
                }
              else
                {
                  insertLocalVectorBatch(cellDim, *workSet, isBCRqc_[r], 
                                         *testLocalDOFs,
                                         group.nTestNodes(),
                                         group.testID(), *localValues, vec);
                }
            }
        }
    }
  mat->freezeValues();

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
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  RefCountPtr<Array<int> > testLocalDOFs 
    = rcp(new Array<int>());

  if (vecNeedsConfiguration_)
    {
      configureVector(b);
    }



  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::assemble()");

  b.zero();

  const Array<EvalContext>& contexts = contexts_.get(VectorOnly);
  const Array<Array<IntegralGroup> >& groups = groups_.get(VectorOnly);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(VectorOnly);

  for (int r=0; r<rqc_.size(); r++)
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
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);


          workSetCounter++;

          mediators_[r]->setCellBatch(workSet, J);

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

          int nTestNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, *testLocalDOFs,
                                       nTestNodes);

          ElementIntegral::invalidateTransformationMatrices();
          for (int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(*J, vectorCoeffs, 
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
                                     group.nTestNodes(),
                                     group.testID(), *localValues, vec);
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

  SUNDANCE_VERB_LOW("Computing functional and gradient"); 

  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  RefCountPtr<Array<int> > testLocalDOFs 
    = rcp(new Array<int>());

  if (vecNeedsConfiguration_)
    {
      configureVector(gradient);
    }



  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(gradient.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::evaluate()");

  gradient.zero();
  double localSum = 0.0;

  const Array<EvalContext>& contexts = contexts_.get(FunctionalAndGradient);
  const Array<Array<IntegralGroup> >& groups 
    = groups_.get(FunctionalAndGradient);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(FunctionalAndGradient);

  for (int r=0; r<rqc_.size(); r++)
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
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);


          workSetCounter++;

          mediators_[r]->setCellBatch(workSet, J);

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

          int nTestNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, *testLocalDOFs,
                                       nTestNodes);


          ElementIntegral::invalidateTransformationMatrices();

          for (int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(*J, vectorCoeffs, 
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
                                         group.nTestNodes(),
                                         group.testID(), *localValues, vec);
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

  SUNDANCE_VERB_LOW(tab << "Computing functional"); 

  SUNDANCE_VERB_MEDIUM(tab << "work set size is " << workSetSize()); 

  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  Array<RefCountPtr<EvalVector> > vectorCoeffs;
  Array<double> constantCoeffs;
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  double localSum = 0.0;

  const Array<EvalContext>& contexts = contexts_.get(FunctionalOnly);
  const Array<Array<IntegralGroup> >& groups 
    = groups_.get(FunctionalOnly);
  const Array<const EvaluatableExpr*>& evalExprs 
    = evalExprs_.get(FunctionalOnly);

  for (int r=0; r<rqc_.size(); r++)
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
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }

          SUNDANCE_VERB_MEDIUM(tab1 << "doing work set " << workSetCounter
                               << " consisting of " 
                               << workSet->size() << " cells");
          SUNDANCE_VERB_EXTREME("cells are " << *workSet);


          workSetCounter++;

          mediators_[r]->setCellBatch(workSet, J);

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

          ElementIntegral::invalidateTransformationMatrices();
          for (int g=0; g<groups[r].size(); g++)
            {
              const IntegralGroup& group = groups[r][g];
              if (!group.evaluate(*J, vectorCoeffs, 
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
                                       const Array<int>& testIndices,
                                       const Array<int>& unkIndices,
                                       int nTestNodes, int nUnkNodes,
                                       const Array<int>& testID, 
                                       const Array<int>& unkID,
                                       const Array<double>& localValues, 
                                       LoadableMatrix<double>* mat) const 
{
  Tabs tab;
  TimeMonitor timer(matInsertTimer());

  SUNDANCE_VERB_HIGH(tab << "inserting local matrix values...");

  static Array<int> skipRow;

  int nCells = workSet.size();
  int rowsPerFunc = nTestNodes*workSet.size();
  int colsPerFunc = nUnkNodes*workSet.size();

  skipRow.resize(rowsPerFunc);

  int highestIndex = lowestRow_ + rowMap_->numLocalDOFs();

  int lowestLocalRow = rowMap_->lowestLocalDOF();

  if (verbosity() > VerbHigh)
    {
      cerr << "isBC " << isBCRqc << endl;
      cerr << "num test nodes = " << nTestNodes << endl;
      cerr << "num unk nodes = " << nUnkNodes << endl;
      cerr << "num rows per function = " << rowsPerFunc << endl;
      cerr << "num cols per function = " << colsPerFunc << endl;
      cerr << "num cells = " << nCells << endl;
      cerr << "testID = " << testID << endl;
      cerr << "unkID = " << unkID << endl;
    }
  for (int t=0; t<testID.size(); t++)
    {
      for (int u=0; u<unkID.size(); u++)
        {
          if (isBCRqc)
            {
              for (int r=0; r<rowsPerFunc; r++)
                {
                  int row = testIndices[r+testID[t]*rowsPerFunc];
                  int localRowIndex = row - lowestLocalRow;
                  skipRow[r] = row < lowestRow_ || row >= highestIndex
                    || !(*isBCRow_)[localRowIndex];
                }
            }
          else
            {
              for (int r=0; r<rowsPerFunc; r++)
                {
                  int row = testIndices[r+testID[t]*rowsPerFunc];
                  int localRowIndex = row - lowestLocalRow;
                  skipRow[r] = row < lowestRow_ || row >= highestIndex
                    || (*isBCRow_)[localRowIndex];
                }
            }
     
          if (verbosity() > VerbHigh)
            {
              cerr << "adding: " << endl;
              int k = 0;
              for (int cell=0; cell<nCells; cell++)
                {
                  for (int r=0; r<nTestNodes; r++)
                    {
                      cerr << testIndices[testID[t]*rowsPerFunc + cell*nTestNodes+r]
                           << ": {";
                      for (int c=0; c<nUnkNodes; c++,k++)
                        {
                          cerr << "(" << unkIndices[unkID[u]*colsPerFunc+cell*nUnkNodes+c]
                               << ", " << localValues[k] << ")";
                          if (c < (nUnkNodes-1)) cerr << ", ";
                        }
                      cerr << "}" << endl;
                    }
                }
            }
          mat->addElementBatch(rowsPerFunc,
                               nTestNodes,
                               &(testIndices[testID[t]*rowsPerFunc]),
                               nUnkNodes,
                               &(unkIndices[unkID[u]*colsPerFunc]),
                               &(localValues[0]),
                               &(skipRow[0]));
        }
    }
}


/* ------------  insert elements into the vector  ------------- */

void Assembler::insertLocalVectorBatch(int cellDim, 
                                       const Array<int>& workSet, 
                                       bool isBCRqc,
                                       const Array<int>& testIndices,
                                       int nTestNodes, 
                                       const Array<int>& testID, 
                                       const Array<double>& localValues, 
                                       TSFExtended::LoadableVector<double>* vec) const 
{
  TimeMonitor timer(vecInsertTimer());
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "inserting local vector values...");
  SUNDANCE_VERB_EXTREME(tab << "values are " << localValues);
  int rowsPerFunc = nTestNodes*workSet.size();

  int lowestLocalRow = rowMap_->lowestLocalDOF();

  for (int i=0; i<testID.size(); i++)
    {
      for (int r=0; r<rowsPerFunc; r++)
        {
          int rowIndex = testIndices[r + testID[i]*rowsPerFunc];
          int localRowIndex = rowIndex - lowestLocalRow;
          if (!(rowMap_->isLocalDOF(rowIndex))
              || isBCRqc!=(*isBCRow_)[localRowIndex]) continue;
          {
            vec->addToElement(rowIndex, localValues[r]);
          }
        }
    }
  SUNDANCE_VERB_HIGH(tab << "...done");
}



/* ------------  get the nonzero pattern for the matrix ------------- */
                       
                       
void Assembler::getGraph(Array<int>& graphData,
                         Array<int>& rowPtrs,
                         Array<int>& nnzPerRow) const 
{
  TimeMonitor timer(graphBuildTimer());
  Tabs tab;




  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());

  RefCountPtr<Array<int> > testLocalDOFs 
    = rcp(new Array<int>());

  RefCountPtr<Array<int> > unkLocalDOFs
    = rcp(new Array<int>());

  SUNDANCE_OUT(verbosity() > VerbLow, tab << "Creating graph: there are " << rowMap()->numLocalDOFs()
               << " local equations");


  int lowestLocalRow = rowMap_->lowestLocalDOF();


  Array<Set<int> > tmpGraph;
  tmpGraph.resize(rowMap()->numLocalDOFs());
  //  int cap = eqn_->numUnks() * (int) round(pow(3.0, mesh_.spatialDim()));
  //  cerr << "capacity = " << cap << endl;
  // {
//     TimeMonitor timer2(tmpGraphBuildTimer());
//     for (int i=0; i<tmpGraph.size(); i++)
//       {
//         tmpGraph[i].setCapacity(cap);
//       }
//   }


  {
    TimeMonitor timer2(colSearchTimer());
    for (int d=0; d<eqn_->numRegions(); d++)
      {
        Tabs tab0;
        CellFilter domain = eqn_->region(d);
        SUNDANCE_OUT(verbosity() > VerbMedium, 
                     tab0 << "cell set " << domain
                     << " isBCRegion=" << eqn_->isBCRegion(d));
        int dim = domain.dimension(mesh_);
        CellSet cells = domain.getCells(mesh_);

        RefCountPtr<Set<OrderedPair<int, int> > > pairs ;
        if (eqn_->hasVarUnkPairs(domain)) pairs = eqn_->varUnkPairs(domain);

        SUNDANCE_OUT(verbosity() > VerbMedium && pairs.get() != 0, 
                     tab0 << "non-BC pairs = "
                     << *pairs);
       
        RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
        if (eqn_->isBCRegion(d))
          {
            if (eqn_->hasBCVarUnkPairs(domain)) 
              {
                bcPairs = eqn_->bcVarUnkPairs(domain);
                SUNDANCE_OUT(verbosity() > VerbMedium, tab0 << "BC pairs = "
                             << *bcPairs);
              }
          }
        Array<Set<int> > unksForTestsSet(eqn_->numVars());
        Array<Set<int> > bcUnksForTestsSet(eqn_->numVars());

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
                bcUnksForTestsSet[eqn_->reducedVarID(t)].put(eqn_->reducedUnkID(u));
              }
          }

        Array<Array<int> > unksForTests(unksForTestsSet.size());
        Array<Array<int> > bcUnksForTests(bcUnksForTestsSet.size());

        for (int t=0; t<unksForTests.size(); t++)
          {
            unksForTests[t] = unksForTestsSet[t].elements();
            bcUnksForTests[t] = bcUnksForTestsSet[t].elements();
          }
      
        int nTestFuncs = 1;
        int nUnkFuncs = 1;
        int nTestNodes;
        int nUnkNodes;

      

        int highestRow = lowestRow_ + rowMap_->numLocalDOFs();

        int owner;
        int nt = eqn_->numVars();
        CellIterator iter=cells.begin();
        while (iter != cells.end())
          {
            /* build a work set */
            workSet->resize(0);
            for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
              {
                workSet->append(*iter);
              }

            int nCells = workSet->size();

            rowMap_->getDOFsForCellBatch(dim, *workSet, *testLocalDOFs,
                                         nTestNodes);
            if (rowMap_.get()==colMap_.get())
              {
                unkLocalDOFs = testLocalDOFs;
                nUnkNodes = nTestNodes;
              }
            else
              {
                colMap_->getDOFsForCellBatch(dim, *workSet, 
                                             *unkLocalDOFs, nUnkNodes);
              }
          
            if (pairs.get() != 0)
              {
                for (int c=0; c<workSet->size(); c++)
                  {
                    int testOffset = c*nTestNodes*nTestFuncs;
                    int unkOffset = c*nUnkNodes*nUnkFuncs;
                    for (int t=0; t<nt; t++)
                      {
                        for (int uit=0; uit<unksForTests[t].size(); uit++)
                          {
                            Tabs tab2;
                            int u = unksForTests[t][uit];
                            for (int n=0; n<nTestNodes; n++)
                              {
                                int row
                                  = (*testLocalDOFs)[(t*nCells + c)*nTestNodes+n];
                                if (row < lowestRow_ || row >= highestRow
                                    || (*isBCRow_)[row-lowestRow_]) continue;
                                Set<int>& colSet = tmpGraph[row-lowestRow_];
                                for (int m=0; m<nUnkNodes; m++)
                                  {
                                    int col 
                                      = (*unkLocalDOFs)[(u*nCells+c)*nUnkNodes+m];
                                    colSet.put(col);
                                  }
                              }
                          }
                      }
                  }
              }
            if (bcPairs.get() != 0)
              {
                for (int c=0; c<workSet->size(); c++)
                  {
                    int testOffset = c*nTestNodes*nTestFuncs;
                    int unkOffset = c*nUnkNodes*nUnkFuncs;
                    for (int t=0; t<nt; t++)
                      {
                        for (int uit=0; uit<bcUnksForTests[t].size(); uit++)
                          {
                            Tabs tab2;
                            int u = bcUnksForTests[t][uit];
                            for (int n=0; n<nTestNodes; n++)
                              {
                                int row
                                  = (*testLocalDOFs)[(t*nCells + c)*nTestNodes+n];
                                if (row < lowestRow_ || row >= highestRow
                                    || !(*isBCRow_)[row-lowestRow_]) continue;
                                Set<int>& colSet = tmpGraph[row-lowestRow_];
                                for (int m=0; m<nUnkNodes; m++)
                                  {
                                    int col 
                                      = (*unkLocalDOFs)[(u*nCells+c)*nUnkNodes+m];
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
    int nLocalRows = rowMap()->numLocalDOFs();

    int nnz = 0;
    rowPtrs.resize(nLocalRows);
    nnzPerRow.resize(rowMap()->numLocalDOFs());
    for (int i=0; i<nLocalRows; i++) 
      {
        rowPtrs[i] = nnz;
        nnzPerRow[i] = tmpGraph[i].size();
        nnz += nnzPerRow[i];
      }

    graphData.resize(nnz);
    int* base = &(graphData[0]);
    for (int i=0; i<nLocalRows; i++)
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
