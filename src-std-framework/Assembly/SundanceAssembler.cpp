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

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


Assembler
::Assembler(const Mesh& mesh, 
            const RefCountPtr<EquationSet>& eqn,
            const VectorType<double>& vectorType,
            const VerbositySetting& verb)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
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
  verbosity() = verb;

  RefCountPtr<GrouperBase> grouper = rcp(new TrivialGrouper());
  grouper->verbosity() = verbosity();

  DOFMapBuilder mapBuilder(mesh, eqn);

  rowMap_ = mapBuilder.rowMap();
  colMap_ = mapBuilder.colMap();
  bcRows_ = mapBuilder.bcRows();

  lowestRow_ = mapBuilder.rowMap()->lowestLocalDOF();

  isBCRow_.resize(colMap()->numLocalDOFs());
  for (int i=0; i<isBCRow_.size(); i++) isBCRow_[i] = false;
  for (Set<int>::const_iterator i=bcRows()->begin(); i != bcRows()->end(); i++)
    {
      int row = *i;
      if (rowMap_->isLocalDOF(row)) isBCRow_[row-lowestRow_] = true;
    } 

  

  for (int r=0; r<eqn->regionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integral groups for rqc=" << rqc);
                         
      rqc_.append(rqc);
      isBCRqc_.append(false);
      const Expr& expr = eqn->expr(rqc);
      const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(rqc);
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_.append(ee);
      const RefCountPtr<SparsityPattern>& sparsity 
        = ee->sparsity(ee->getDerivSetIndex(rqc));
      Array<IntegralGroup> groups;
      grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
      groups_.append(groups);
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }

  
  for (int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integral groups for BC rqc=" << rqc);
      rqc_.append(rqc);
      isBCRqc_.append(true);
      const Expr& expr = eqn->bcExpr(rqc);
      const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(rqc);
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());
      const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
      evalExprs_.append(ee);
      const RefCountPtr<SparsityPattern>& sparsity 
        = ee->sparsity(ee->getDerivSetIndex(rqc));
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

  VectorSpace<double> rowSpace = vecType_.createSpace(rowMap()->numDOFs(),
                                                      rowMap()->numLocalDOFs(),
                                                      &(localRowIndices[0]));
  VectorSpace<double> colSpace = vecType_.createSpace(colMap()->numDOFs(),
                                                      colMap()->numLocalDOFs(),
                                                      &(localColIndices[0]));

  
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
  Array<Set<int> > graph;
  Array<int> colIndices;
  Array<double> zeros;
  getGraph(graph);

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "...done");
  if (verbosity() > VerbHigh)
    {
      Tabs tab1;
      cerr << tab1 << "graph: " << endl;
      for (int i=0; i<graph.size(); i++) 
        {
          Tabs tab2;
          cerr << tab2 << "row=" << i << " " << graph[i] << endl;
        }
    }

  
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: initializing matrix and vector...");
  for (int i=0; i<graph.size(); i++)
    {
      colIndices.resize(graph[i].size());
      zeros.resize(graph[i].size());
      Set<int>::const_iterator iter;
      int j=0;
      for (iter=graph[i].begin(); iter != graph[i].end(); iter++, j++)
        {
          colIndices[j] = *iter;
          zeros[j] = 0.0;
        }
      mat->setRowValues(lowestRow_ + i, colIndices.size(),
                        &(colIndices[0]), &(zeros[0]));
      vec->setElement(lowestRow_ + i, 0.0);
    }
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

  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());
  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  RefCountPtr<EvalVectorArray> coeffs;
  RefCountPtr<CellJacobianBatch> J = rcp(new CellJacobianBatch());

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

          for (int g=0; g<groups_[r].size(); g++)
            {
              const IntegralGroup& group = groups_[r][g];
              if (!group.evaluate(*J, coeffs, localValues)) continue;
              
              if (group.isTwoForm())
                {
                  insertLocalMatrixValues(cellDim, *workSet, isBCRqc_[r],
                                          group.nTestNodes(),
                                          group.nUnkNodes(),
                                          group.testID(), group.unkID(), 
                                          *localValues, mat);
                }
              else
                {
                  insertLocalVectorValues(cellDim, *workSet, isBCRqc_[r], 
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

  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: inserting local matrix values");
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab << "Assembler: values are " << localValues);

  Array<int> testIndices;
  Array<int> unkIndices;
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
          Array<double> x(unkIndices.size());
          Array<int> p(unkIndices.size());
          for (int r=0; r<testIndices.size(); r++)
            {
              if (rowMap_->isLocalDOF(testIndices[r]) && 
                  isBCRqc==isBCRow_[testIndices[r]])
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

void Assembler::insertLocalVectorValues(int cellDim, 
                                        const Array<int>& workSet, 
                                        bool isBCRqc,
                                        int nTestNodes, 
                                        const Array<int>& testID, 
                                        const Array<double>& localValues, 
                                        TSFExtended::LoadableVector<double>* vec) const 
{

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
                  isBCRqc==isBCRow_[testIndices[r]])
                {
                  vec->addToElement(testIndices[r], localValues[c*nTestNodes+r]);
                }
            }
        }
    }
}

                       
                       
void Assembler::getGraph(Array<Set<int> >& graph) const 
{
  Tabs tab;

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
      Array<Set<int> > unksForTests(eqn_->numTests());
      Array<Set<int> > bcUnksForTests(eqn_->numTests());

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

              unksForTests[eqn_->reducedTestID(t)].put(eqn_->reducedUnkID(u));
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
              bcUnksForTests[eqn_->reducedTestID(t)].put(eqn_->reducedUnkID(u));
            }
        }
      
      Array<int> testDOFs;
      Array<int> unkDOFs;
      int owner;
      int nt = eqn_->numTests();
      for (CellIterator cell=cells.begin(); cell != cells.end(); cell++)
        {
          int cellLID = *cell;
          Tabs tab1;
          SUNDANCE_OUT(verbosity() > VerbMedium, tab1 
                       << "cell dim=" << dim << " cell=" << cellLID);
          if (pairs.get() != 0)
            {
              for (int t=0; t<nt; t++)
                {
                  if (unksForTests[t].size()==0) continue;
                  rowMap()->getDOFsForCell(dim, cellLID, t, testDOFs);
                  for (Set<int>::const_iterator uit=unksForTests[t].begin();
                       uit != unksForTests[t].end(); uit++)
                    {
                      Tabs tab2;
                      int u = *uit;
                      colMap()->getDOFsForCell(dim, cellLID, u, unkDOFs);
                      for (int r=0; r<testDOFs.length(); r++)
                        {
                          int row = testDOFs[r];
                          SUNDANCE_OUT(verbosity() > VerbMedium,
                                       tab2 << "row " << row 
                                       << " isBC=" << isBCRow(row)
                                       << " isLocal=" 
                                       << colMap()->isLocalDOF(row));
                          if (!colMap()->isLocalDOF(row) || isBCRow(row)) continue;
                          for (int c=0; c<unkDOFs.length(); c++)
                            {
                              graph[row-lowestRow_].put(unkDOFs[c]);
                            }
                        }
                    }
                }
            }
          if (bcPairs.get() != 0)
            {
              for (int t=0; t<nt; t++)
                {
                  if (bcUnksForTests[t].size()==0) continue;
                  rowMap()->getDOFsForCell(dim, cellLID, t, testDOFs);
                  for (Set<int>::const_iterator uit=bcUnksForTests[t].begin();
                       uit != bcUnksForTests[t].end(); uit++)
                    {
                      Tabs tab2;
                      int u = *uit;
                      colMap()->getDOFsForCell(dim, cellLID, u, unkDOFs);
                      for (int r=0; r<testDOFs.length(); r++)
                        {
                          int row = testDOFs[r];
                          SUNDANCE_OUT(verbosity() > VerbMedium,
                                       tab2 << "row " << row 
                                       << " isBC=" << isBCRow(row)
                                       << " isLocal=" 
                                       << colMap()->isLocalDOF(row));
                          if (!colMap()->isLocalDOF(row) || !isBCRow(row)) continue;
                          for (int c=0; c<unkDOFs.length(); c++)
                            {
                              graph[row-lowestRow_].put(unkDOFs[c]);
                            }
                        }
                    }
                }
            }
        }
    }
}
