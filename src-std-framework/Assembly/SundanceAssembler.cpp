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
    contexts_(2),
    isBCRqc_(),
    groups_(2),
    mediators_(),
    evalExprs_(2),
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

      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());

      for (int order=1; order<=2; order++)
        {
          const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(order, rqc);
          EvalContext context = eqn->rqcToContext(order, rqc);
          contexts_[order-1].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[order-1].append(ee);
          const RefCountPtr<SparsityPattern>& sparsity 
            = ee->sparsity(ee->getDerivSetIndex(context));
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       "sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
          groups_[order-1].append(groups);
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

      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integral groups for rqc=" << rqc << endl
                   << "expr = " << expr.toXML().toString());
      
      int cellDim = CellFilter(rqc.domain()).dimension(mesh);
      CellType cellType = mesh.cellType(cellDim);
      QuadratureFamily quad(rqc.quad());

      for (int order=1; order<=2; order++)
        {
          const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(order, rqc);
          EvalContext context = eqn->bcRqcToContext(order, rqc);
          contexts_[order-1].append(context);
          const EvaluatableExpr* ee = EvaluatableExpr::getEvalExpr(expr);
          evalExprs_[order-1].append(ee);
          const RefCountPtr<SparsityPattern>& sparsity 
            = ee->sparsity(ee->getDerivSetIndex(context));
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       "sparsity pattern " << *sparsity);

          Array<IntegralGroup> groups;
          grouper->findGroups(*eqn, cellType, cellDim, quad, sparsity, groups);
          groups_[order-1].append(groups);
        }
      mediators_.append(rcp(new QuadratureEvalMediator(mesh, cellDim, 
                                                       quad)));
    }
}

void Assembler::configureVector(Vector<double>& b) const 
{
  Tabs tab;
  TimeMonitor timer(configTimer());
  VectorSpace<double> colSpace = colSpace_->vecSpace();
  
  b = colSpace.createMember();

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

  A = vecType_.createMatrix(colSpace, rowSpace);

  if (vecNeedsConfiguration_)
    {
      configureVector(b);
    }
                                                      
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

  b.zero();
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "...done");

  

  matNeedsConfiguration_ = false;
}


/* ------------  assemble both the vector and the matrix  ------------- */

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
      cerr << "map" << endl;
      rowMap_->print(cerr);
      cerr << "BC row flags " << endl;
      for (int i=0; i<isBCRow_->size(); i++) 
        {
          cerr << i << " " << (*isBCRow_)[i] << endl;
        }
    }

  for (int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tab0 << "doing matrix/vector assembly for rqc=" << rqc_[r]);

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts_[1][r]);

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

          evalExprs_[1][r]->flushResultCache();
          mesh_.getJacobians(cellDim, *workSet, *J);
          evalExprs_[1][r]->evaluate(*evalMgr_, coeffs);

          int nTestNodes;
          int nUnkNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, *testLocalDOFs,
                                       nTestNodes);
          SUNDANCE_OUT(verbosity() > VerbHigh,
                       tab1 << "local DOF values " << *testLocalDOFs);
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
          
          for (int g=0; g<groups_[1][r].size(); g++)
            {
              const IntegralGroup& group = groups_[1][r][g];
              if (!group.evaluate(*J, coeffs, localValues)) continue;

              if (verbosity() > VerbHigh)
                {
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

  if (verbosity() > VerbHigh)
    {
      cerr << "matrix = " << endl;
      A.print(cerr);
    }

}


/* ------------  assemble the vector alone  ------------- */

void Assembler::assemble(Vector<double>& b) const 
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

  if (vecNeedsConfiguration_)
    {
      configureVector(b);
    }



  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(b.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in Assembler::assemble()");

  b.zero();

  for (int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tab0 << "doing matrix/vector assembly for rqc=" << rqc_[r]);

      /* specify the mediator for this RQC */
      evalMgr_->setMediator(mediators_[r]);
      evalMgr_->setRegion(contexts_[0][r]);

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

          evalExprs_[0][r]->flushResultCache();
          mesh_.getJacobians(cellDim, *workSet, *J);
          evalExprs_[0][r]->evaluate(*evalMgr_, coeffs);

          int nTestNodes;
          rowMap_->getDOFsForCellBatch(cellDim, *workSet, *testLocalDOFs,
                                       nTestNodes);


          for (int g=0; g<groups_[0][r].size(); g++)
            {
              const IntegralGroup& group = groups_[0][r][g];
              if (!group.evaluate(*J, coeffs, localValues)) continue;
              
              insertLocalVectorBatch(cellDim, *workSet, isBCRqc_[r], 
                                     *testLocalDOFs,
                                     group.nTestNodes(),
                                     group.testID(), *localValues, vec);
            }
        }
    }
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

  static Array<int> skipRow;

  int rowsPerFunc = nTestNodes*workSet.size();
  int colsPerFunc = nUnkNodes*workSet.size();

  skipRow.resize(rowsPerFunc);

  int highestIndex = lowestRow_ + rowMap_->numLocalDOFs();

  if (verbosity() > VerbHigh)
    {
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
                  skipRow[r] = row < lowestRow_ || row >= highestIndex
                    || !(*isBCRow_)[row];
                }
            }
          else
            {
              for (int r=0; r<rowsPerFunc; r++)
                {
                  int row = testIndices[r+testID[t]*rowsPerFunc];
                  skipRow[r] = row < lowestRow_ || row >= highestIndex
                    || (*isBCRow_)[row];
                }
            }
     
          if (verbosity() > VerbHigh)
            {
              cerr << "adding: " << endl;
              int k = 0;
              for (int r=0; r<nTestNodes; r++)
                {
                  cerr << testIndices[testID[t]*rowsPerFunc+r]
                       << ": {";
                  for (int c=0; c<nUnkNodes; c++,k++)
                    {
                      cerr << "(" << unkIndices[unkID[u]*colsPerFunc+c]
                           << ", " << localValues[k] << ")";
                      if (c < (nUnkNodes-1)) cerr << ", ";
                    }
                  cerr << "}" << endl;
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
  SUNDANCE_OUT(verbosity() > VerbLow, 
               tab << "Assembler: inserting local vector values");
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               tab << "Assembler: values are " << localValues);
  int rowsPerFunc = nTestNodes*workSet.size();

  for (int i=0; i<testID.size(); i++)
    {
      for (int r=0; r<rowsPerFunc; r++)
        {
          int rowIndex = testIndices[r + testID[i]*rowsPerFunc];
          if (isBCRqc!=(*isBCRow_)[rowIndex] 
              || !(rowMap_->isLocalDOF(rowIndex))) continue;
          {
            vec->addToElement(rowIndex, localValues[r]);
          }
        }
    }
}



/* ------------  get the nonzero pattern for the matrix ------------- */
                       
                       
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
                          const int* rowPtr 
                            = &((*testLocalDOFs)[testOffset]);
                          for (int n=0; n<nTestNodes; n++)
                            {
                              int row = rowPtr[t*nTestNodes+n];
                              // if (!rowMap()->isLocalDOF(row) 
//                                   || isBCRow(row)) continue;
                              if (row < lowestRow_ || row >= highestRow
                                  || (*isBCRow_)[row]) continue;
                              ColSetType& colSet = graph[row-lowestRow_];
                              const int* colPtr 
                                = &((*unkLocalDOFs)[unkOffset]);
                              for (int m=0; m<nUnkNodes; m++)
                                {
                                  int col = colPtr[u*nUnkNodes+m];
                                  colSet.insert(col);
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
                              int row = (*testLocalDOFs)[testOffset + t*nTestNodes + n];
                              // if (!rowMap()->isLocalDOF(row) 
//                                   || !isBCRow(row)) continue;
                              if (row < lowestRow_ || row >= highestRow
                                  || !(*isBCRow_)[row]) continue;
                              ColSetType& colSet = graph[row-lowestRow_];
                              for (int m=0; m<nUnkNodes; m++)
                                {
                                  int col = (*unkLocalDOFs)[unkOffset + u*nUnkNodes + m];
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
