/* @HEADER@ */
/* @HEADER@ */

#include "SundanceAssembler.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceTrivialGrouper.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

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


Assembler
::Assembler(const Mesh& mesh, 
            const RefCountPtr<EquationSet>& eqn,
            const VectorType<double>& vectorType,
            const VerbositySetting& verb)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    rowSpace_(),
    colSpace_(),
    bcRows_(),
    rqc_(),
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

  RefCountPtr<GrouperBase> grouper = rcp(new TrivialGrouper());
  grouper->verbosity() = verbosity();

  DOFMapBuilder mapBuilder(mesh, eqn);

  rowMap_ = mapBuilder.rowMap();
  colMap_ = mapBuilder.colMap();
  rowSpace_ = rcp(new DiscreteSpace(mesh, mapBuilder.testBasisArray(), rowMap_, vecType_));
  colSpace_ = rcp(new DiscreteSpace(mesh, mapBuilder.unkBasisArray(), colMap_, vecType_));
  isBCRow_ = mapBuilder.isBCRow();

  lowestRow_ = mapBuilder.rowMap()->lowestLocalDOF();


  for (int r=0; r<eqn->regionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];
                         
      rqc_.append(rqc);
      isBCRqc_.append(false);
      const Expr& expr = eqn->expr(rqc);

      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integral groups for rqc=" << rqc << endl
                   << "expr = " << expr);
      const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(rqc);
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_.append(ee);
      const RefCountPtr<SparsityPattern>& sparsity 
        = ee->sparsity(ee->getDerivSetIndex(rqc));
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "sparsity pattern " << *sparsity);
      
      Array<IntegralGroup> groups;
      grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
      groups_.append(groups);
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }

  
  for (int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];
      rqc_.append(rqc);
      isBCRqc_.append(true);
      const Expr& expr = eqn->bcExpr(rqc);

      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integral groups for rqc=" << rqc << endl
                   << "expr = " << expr.toXML().toString());
      const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(rqc);
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_.append(ee);
      const RefCountPtr<SparsityPattern>& sparsity 
        = ee->sparsity(ee->getDerivSetIndex(rqc));
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "sparsity pattern " << *sparsity);
      Array<IntegralGroup> groups;
      grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
      groups_.append(groups);
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }
}


void Assembler::configureMat(LinearOperator<double>& A,
                             Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());
  
  Array<int> localRowIndices(rowMap()->numLocalDOFs());
  for (int i=0; i<localRowIndices.size(); i++) 
    {
      localRowIndices[i] = lowestRow_+i;
    }
  Array<int> localColIndices(colMap()->numLocalDOFs());
  int lowestCol = colMap()->lowestLocalDOF();
  for (int i=0; i<localColIndices.size(); i++) 
    {
      localColIndices[i] = lowestCol+i;
    }

  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab << "Assembler: num rows = " << rowMap()->numDOFs());
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               tab << "Assembler: num cols = " << colMap()->numDOFs());

  
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: creating row and col spaces...");

  VectorSpace<double> rowSpace = rowSpace_->vecSpace();
  VectorSpace<double> colSpace = colSpace_->vecSpace();

  
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "...done");

  b = colSpace.createMember();
  A = vecType_.createMatrix(colSpace, rowSpace);
                                                      
  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::assemble()");

  TSFExtended::LoadableMatrix<double>* mat
    = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(A.ptr().get());

  TEST_FOR_EXCEPTION(mat==0, RuntimeError,
                     "matrix is not loadable in Assembler::assemble()");


  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: creating graph...");
  Array<ColSetType<int> > graph;
  Array<int> colIndices;
  Array<double> zeros;
  getGraph(graph);

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "...done");
  

  
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: initializing matrix and vector...");

  mat->configure(lowestRow_, graph);

  // int maxSize = 0;
//   for (int i=0; i<graph.size(); i++)
//     {
//       graph[i].elements(colIndices);
//       if (colIndices.size() > maxSize) 
//         {
//           zeros.resize(colIndices.size());
//           for (int j=maxSize; j<zeros.size(); j++) zeros[j] = 0.0;
//           maxSize = zeros.size();
//         }
//       mat->setRowValues(lowestRow_ + i, colIndices.size(),
//                         &(colIndices[0]), &(zeros[0]));
//     }

  b.zero();
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "...done");

  if (verbosity() > VerbHigh)
    {
      Tabs tab1;
      cerr << tab1 << "Assemble: matrix before fill = ";
      A.print(cerr);
      cerr << tab1 << "Assemble: vector before fill = ";
      A.print(cerr);
    }
}


void Assembler::assemble(LinearOperator<double>& A,
                         Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(assemblyTimer());

  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());
  SUNDANCE_OUT(verbosity() > VerbSilent, 
               "work set size is " << workSetSize()); 
  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  RefCountPtr<EvalVectorArray> coeffs;
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

  RefCountPtr<Array<int> > testLocalDOFs 
    = rcp(new Array<int>());

  RefCountPtr<Array<int> > unkLocalDOFs
    = rcp(new Array<int>());

  configureMat(A, b);

  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::assemble()");

  TSFExtended::LoadableMatrix<double>* mat
    = dynamic_cast<TSFExtended::LoadableMatrix<double>* >(A.ptr().get());

  TEST_FOR_EXCEPTION(mat==0, RuntimeError,
                     "matrix is not loadable in Assembler::assemble()");

  for (int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tab0 << "doing matrix/vector assembly for rqc=" << rqc_[r]);

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(rqc_[r]);

      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);
      mediators_[r]->setCellType(cellType);      

      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tab0 << "cell type = " << cellType);

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
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tab1 << "doing work set " << workSetCounter
                       << " consisting of " << workSet->size() << " cells");
          workSetCounter++;

          mediators_[r]->setCellBatch(workSet);

          evalExprs_[r]->flushResultCache();
          mesh_.getJacobians(cellDim, *workSet, *J);
          evalExprs_[r]->evaluate(*evalMgr_, coeffs);

          int nTestNodes;
          int nUnkNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, 0, *testLocalDOFs,
                                       nTestNodes);
          if (rowMap_.get()==colMap_.get())
            {
              unkLocalDOFs = testLocalDOFs;
              nUnkNodes = nTestNodes;
            }
          else
            {
              colMap_->getDOFsForCellBatch(cellDim, *workSet, 0, 
                                           *unkLocalDOFs,
                                           nUnkNodes);
            }

          for (int g=0; g<groups_[r].size(); g++)
            {
              const IntegralGroup& group = groups_[r][g];
              if (!group.evaluate(*J, coeffs, localValues)) continue;
              
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

}

void Assembler::insertLocalMatrixValues(int cellDim, 
                                        const Array<int>& workSet, 
                                        bool isBCRqc,
                                        int nTestNodes, int nUnkNodes,
                                        const Array<int>& testID, 
                                        const Array<int>& unkID,
                                        const Array<double>& localValues, 
                                        LoadableMatrix<double>* mat) const 
{
  Tabs tab;
  TimeMonitor timer(matInsertTimer());

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: inserting local matrix values");
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab << "Assembler: values are " << localValues);

  Array<int> testIndices;
  Array<int> unkIndices;
  Array<double> x(nUnkNodes);
  Array<int> p(nUnkNodes);
  int nNodes = nTestNodes*nUnkNodes;

  SUNDANCE_OUT(verbosity() > VerbMedium, 
               tab << "Assembler: num nodes test=" << nTestNodes
               << " unk=" << nUnkNodes);

  for (int c=0; c<workSet.size(); c++)
    {
      Tabs tab1;
      SUNDANCE_OUT(verbosity() > VerbHigh, 
                   tab1 << "cell=" << c);
      for (int i=0; i<testID.size(); i++)
        { 
          Tabs tab2;
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       tab2 << "testID=" << testID[i]);
          rowMap_->getDOFsForCell(cellDim, workSet[c], testID[i], testIndices);
          colMap_->getDOFsForCell(cellDim, workSet[c], unkID[i], unkIndices);
          
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       tab2 << "row indices=" << testIndices);
          
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       tab2 << "col indices=" << unkIndices);
          //          Array<double> x(unkIndices.size());
          //          Array<int> p(unkIndices.size());
          for (int r=0; r<testIndices.size(); r++)
            {
              if (rowMap_->isLocalDOF(testIndices[r]) && 
                  isBCRqc==(*isBCRow_)[testIndices[r]])
                {
                  for (int j=0; j<x.size(); j++)
                    {
                      p[j] = c*nNodes + j*nUnkNodes + r;
                      x[j] = localValues[p[j]];
                    }
                  const double* data = &(localValues[c*nNodes + r*nUnkNodes]);
                  mat->addToRow(testIndices[r], unkIndices.size(), 
                                &(unkIndices[0]), &(x[0]));
                }
            }
        }
    }
}

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

  static Array<int> skipRow;

  skipRow.resize(testIndices.size());

  int highestIndex = lowestRow_ + rowMap_->numLocalDOFs();

  if (isBCRqc)
    {
      for (int r=0; r<testIndices.size(); r++)
        {
          int row = testIndices[r];
          skipRow[r] = row < lowestRow_ || row >= highestIndex
            || !(*isBCRow_)[row];
        }
    }
  else
    {
      for (int r=0; r<testIndices.size(); r++)
        {
          int row = testIndices[r];
          skipRow[r] = row < lowestRow_ || row >= highestIndex
            || (*isBCRow_)[row];
        }
    }

  mat->addElementBatch(testIndices.size(),
                       nTestNodes,
                       &(testIndices[0]),
                       nUnkNodes,
                       &(unkIndices[0]),
                       &(localValues[0]),
                       &(skipRow[0]));

  // SUNDANCE_OUT(verbosity() > VerbLow, 
//                tab << "Assembler: inserting local matrix values");
//   SUNDANCE_OUT(verbosity() > VerbHigh, 
//                tab << "Assembler: values are " << localValues);

//   Array<double> x(nUnkNodes);
//   Array<int> p(nUnkNodes);
//   Array<int> col(nUnkNodes);
//   int nNodes = nTestNodes*nUnkNodes;

//   SUNDANCE_OUT(verbosity() > VerbMedium, 
//                tab << "Assembler: num nodes test=" << nTestNodes
//                << " unk=" << nUnkNodes);

//   for (int i=0; i<testID.size(); i++)
//     {
//       for (int c=0; c<workSet.size(); c++)
//         {
//           for (int r=0; r<nTestNodes; r++)
//             {
//               int rowIndex = testIndices[r + nTestNodes*c];
//               if (isBCRqc!=(*isBCRow_)[rowIndex] 
//                   || !(rowMap_->isLocalDOF(rowIndex))) continue;
//               for (int j=0; j<nUnkNodes; j++)
//                 {
//                   col[j] = unkIndices[j + nUnkNodes*c];
//                   x[j] = localValues[c*nNodes + j*nUnkNodes + r];
//                 }
//               mat->addToRow(rowIndex, nUnkNodes,
//                             &(col[0]), &(x[0]));
//             }
//         }
//     }
}


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
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: inserting local matrix values");
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab << "Assembler: values are " << localValues);

  for (int i=0; i<testID.size(); i++)
    {
      for (int c=0; c<workSet.size(); c++)
        {
          for (int r=0; r<nTestNodes; r++)
            {
              int rowIndex = testIndices[r + nTestNodes*c];
              if (isBCRqc!=(*isBCRow_)[rowIndex] 
                  || !(rowMap_->isLocalDOF(rowIndex))) continue;
              {
                vec->addToElement(rowIndex, localValues[c*nTestNodes+r]);
              }
            }
        }
    }
}


void Assembler::insertLocalVectorValues(int cellDim, 
                                        const Array<int>& workSet, 
                                        bool isBCRqc,
                                        int nTestNodes, 
                                        const Array<int>& testID, 
                                        const Array<double>& localValues, 
                                        TSFExtended::LoadableVector<double>* vec) const 
{
  TimeMonitor timer(vecInsertTimer());
  Tabs tab;
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: inserting local matrix values");
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab << "Assembler: values are " << localValues);

  Array<int> testIndices;
  
  for (int c=0; c<workSet.size(); c++)
    {
      Tabs tab1;
      SUNDANCE_OUT(verbosity() > VerbHigh, 
                   tab1 << "cell=" << c);
      for (int i=0; i<testID.size(); i++)
        {
          Tabs tab2;
          rowMap_->getDOFsForCell(cellDim, workSet[c], testID[i], testIndices);

          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       tab2 << "row indices=" << testIndices);

          for (int r=0; r<testIndices.size(); r++)
            {
              if (rowMap_->isLocalDOF(testIndices[r]) && 
                  isBCRqc==(*isBCRow_)[testIndices[r]])
                {
                  vec->addToElement(testIndices[r], localValues[c*nTestNodes+r]);
                }
            }
        }
    }
}

                       
                       
void Assembler::getGraph(Array<ColSetType>& graph) const 
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

  graph.resize(rowMap()->numLocalDOFs());

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
      if (eqn_->hasTestUnkPairs(domain)) pairs = eqn_->testUnkPairs(domain);

      SUNDANCE_OUT(verbosity() > VerbMedium && pairs.get() != 0, 
                   tab0 << "non-BC pairs = "
                   << *pairs);
       
      RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
      if (eqn_->isBCRegion(d))
        {
          if (eqn_->hasBCTestUnkPairs(domain)) 
            {
              bcPairs = eqn_->bcTestUnkPairs(domain);
              SUNDANCE_OUT(verbosity() > VerbMedium, tab0 << "BC pairs = "
                           << *bcPairs);
            }
        }
      Array<Set<int> > unksForTestsSet(eqn_->numTests());
      Array<Set<int> > bcUnksForTestsSet(eqn_->numTests());

      Set<OrderedPair<int, int> >::const_iterator i;
      
      if (pairs.get() != 0)
        {
          for (i=pairs->begin(); i!=pairs->end(); i++)
            {
              const OrderedPair<int, int>& p = *i;
              int t = p.first();
              int u = p.second();

              TEST_FOR_EXCEPTION(!eqn_->hasTestID(t), InternalError,
                                 "Test function ID " << t << " does not appear "
                                 "in equation set");
              TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
                                 "Unk function ID " << u << " does not appear "
                                 "in equation set");

              unksForTestsSet[eqn_->reducedTestID(t)].put(eqn_->reducedUnkID(u));
            }
        }
      if (bcPairs.get() != 0)
        {
          for (i=bcPairs->begin(); i!=bcPairs->end(); i++)
            {
              const OrderedPair<int, int>& p = *i;
              int t = p.first();
              int u = p.second();
              TEST_FOR_EXCEPTION(!eqn_->hasTestID(t), InternalError,
                                 "Test function ID " << t << " does not appear "
                                 "in equation set");
              TEST_FOR_EXCEPTION(!eqn_->hasUnkID(u), InternalError,
                                 "Unk function ID " << u << " does not appear "
                                 "in equation set");
              bcUnksForTestsSet[eqn_->reducedTestID(t)].put(eqn_->reducedUnkID(u));
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
      int nt = eqn_->numTests();
      CellIterator iter=cells.begin();
      while (iter != cells.end())
        {
          /* build a work set */
          workSet->resize(0);
          for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
            {
              workSet->append(*iter);
            }


          rowMap_->getDOFsForCellBatch(dim, *workSet, 0, *testLocalDOFs,
                                       nTestNodes);
          if (rowMap_.get()==colMap_.get())
            {
              unkLocalDOFs = testLocalDOFs;
              nUnkNodes = nTestNodes;
            }
          else
            {
              colMap_->getDOFsForCellBatch(dim, *workSet, 0, 
                                           *unkLocalDOFs, nUnkNodes);
            }
          
          if (pairs.get() != 0)
            {
              for (int t=0; t<nt; t++)
                {
                  for (int uit=0; uit<unksForTests[t].size(); uit++)
                    {
                      Tabs tab2;
                      int u = unksForTests[t][uit];
                      for (int c=0; c<workSet->size(); c++)
                        {
                          int testOffset = c*nTestNodes*nTestFuncs;
                          int unkOffset = c*nUnkNodes*nUnkFuncs;
                          for (int n=0; n<nTestNodes; n++)
                            {
                              int row = (*testLocalDOFs)[testOffset + n*nTestFuncs + t];
                              // if (!rowMap()->isLocalDOF(row) 
//                                   || isBCRow(row)) continue;
                              if (row < lowestRow_ || row >= highestRow
                                  || (*isBCRow_)[row]) continue;
                              std::set<int>& colSet = graph[row-lowestRow_];
                              for (int m=0; m<nUnkNodes; m++)
                                {
                                  int col = (*unkLocalDOFs)[unkOffset + m*nUnkFuncs + t];
                                  colSet.insert(col);
                                }
                            }
                        }
                    }
                }
            }
          if (bcPairs.get() != 0)
            {
              for (int t=0; t<nt; t++)
                {
                  for (int uit=0; uit<bcUnksForTests[t].size(); uit++)
                    {
                      Tabs tab2;
                      int u = bcUnksForTests[t][uit];
                      for (int c=0; c<workSet->size(); c++)
                        {
                          int testOffset = c*nTestNodes*nTestFuncs;
                          int unkOffset = c*nUnkNodes*nUnkFuncs;
                          for (int n=0; n<nTestNodes; n++)
                            {
                              int row = (*testLocalDOFs)[testOffset + n*nTestFuncs + t];
                              // if (!rowMap()->isLocalDOF(row) 
//                                   || !isBCRow(row)) continue;
                              if (row < lowestRow_ || row >= highestRow
                                  || !(*isBCRow_)[row]) continue;
                              std::set<int>& colSet = graph[row-lowestRow_];
                              for (int m=0; m<nUnkNodes; m++)
                                {
                                  int col = (*unkLocalDOFs)[unkOffset + m*nUnkFuncs + t];
                                  colSet.insert(col);
                                  
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
