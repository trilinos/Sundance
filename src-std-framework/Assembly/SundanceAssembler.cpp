/* @HEADER@ */
/* @HEADER@ */

#include "SundanceAssembler.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"

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
            const RefCountPtr<InserterFactoryBase>& inserterFactory,
            const RefCountPtr<IntegratorFactoryBase>& integratorFactory,
            const VectorType<double>& vectorType,
            const VerbositySetting& verb)
  : mesh_(mesh),
    eqn_(eqn),
    rqc_(),
    isBCRqc_(),
    inserter_(),
    integrator_(),
    evalMgr_(rcp(new EvalManager())),
    isBCRow_(),
    lowestRow_()
{
  verbosity() = verb;

  DOFMapBuilder mapBuilder(mesh, eqn);

  inserter_ = inserterFactory->create(mapBuilder.rowMap(), 
                                      mapBuilder.colMap(),
                                      mapBuilder.bcRows());
  lowestRow_ = mapBuilder.rowMap()->lowestLocalDOF();

  isBCRow_.resize(colMap()->numLocalDOFs());
  for (int i=0; i<isBCRow_.size(); i++) isBCRow_[i] = false;
  for (Set<int>::const_iterator i=bcRows()->begin(); i != bcRows()->end(); i++)
    {
      int row = *i;
      if (colMap()->isLocalDOF(row)) isBCRow_[row-lowestRow_] = true;
    } 

  

  for (int r=0; r<eqn->regionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->regionQuadCombos()[r];
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integrator for rqc=" << rqc);
                         
      rqc_.append(rqc);
      isBCRqc_.append(false);
      const Expr& expr = eqn->expr(rqc);
      const DerivSet& derivs = eqn->nonzeroFunctionalDerivs(rqc);
      RefCountPtr<IntegratorBase> integrator 
        = integratorFactory->create(mesh, expr, rqc, derivs, evalMgr_);
      integrator_.append(integrator);
    }

  
  for (int r=0; r<eqn->bcRegionQuadCombos().size(); r++)
    {
      const RegionQuadCombo& rqc = eqn->bcRegionQuadCombos()[r];
      SUNDANCE_OUT(verbosity() > VerbMedium,
                   "creating integrator for BC rqc=" << rqc);
      rqc_.append(rqc);
      isBCRqc_.append(true);
      const Expr& expr = eqn->bcExpr(rqc);
      const DerivSet& derivs = eqn->nonzeroBCFunctionalDerivs(rqc);
      RefCountPtr<IntegratorBase> integrator 
        = integratorFactory->create(mesh, expr, rqc, derivs, evalMgr_);
      integrator_.append(integrator);
    }
}




void Assembler::assemble(LinearOperator<double>& A,
                         Vector<double>& b) const 
{
  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());
  RefCountPtr<LocalMatrixContainer> localMat ;

  for (int r=0; r<rqc_.size(); r++)
    {
      Tabs tab0;
      SUNDANCE_OUT(verbosity() > VerbLow, 
                   tab0 << "doing matrix assembly for rqc=" << rqc_[r]);
      
      /* get the cells for the current domain */
      CellFilter filter = rqc_[r].domain();
      CellSet cells = filter.getCells(mesh_);
      int cellDim = filter.dimension(mesh_);
      CellType cellType = mesh_.cellType(cellDim);

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
                       tab1 << "doing work set " << workSetCounter);
          workSetCounter++;
          /* integrate the work set to obtain local matrices */
          integrator_[r]->integrate(workSet, localMat);
          /* insert the local matrices in the global matrix */
          inserter_->insert(cellDim, workSet, isBCRqc_[r], localMat.get(),
                            A, b);
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
      RefCountPtr<Set<OrderedPair<int, int> > > pairs;
      
      pairs = eqn_->testUnkPairs(domain);
      SUNDANCE_OUT(verbosity() > VerbMedium && pairs.get() != 0, 
                   tab0 << "non-BC pairs = "
                   << *pairs);
       
      RefCountPtr<Set<OrderedPair<int, int> > > bcPairs ;
      if (eqn_->isBCRegion(d))
        {
          bcPairs = eqn_->bcTestUnkPairs(domain);
          SUNDANCE_OUT(verbosity() > VerbMedium, tab0 << "BC pairs = "
                       << *bcPairs);
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
