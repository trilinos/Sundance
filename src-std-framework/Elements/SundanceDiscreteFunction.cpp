/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFunction.hpp"
#include "SundanceHomogeneousDOFMap.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector())
{}


DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Expr& expr,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector())
{
  Tabs tab;

  RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
  workSet->reserve(workSetSize());
  RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

  TSFExtended::LoadableVector<double>* vec 
    = dynamic_cast<TSFExtended::LoadableVector<double>* >(vector_.ptr().get());

  TEST_FOR_EXCEPTION(vec==0, RuntimeError,
                     "vector is not loadable in DiscreteFunction ctor");


  CellFilter filter = new MaximalCellFilter();

  CellSet cells = filter.getCells(space_.mesh());
  int cellDim = filter.dimension(mesh_);
  CellType cellType = space_.mesh().cellType(cellDim);

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

      mediators->setCellBatch(workSet);

      evalExpr->flushResultCache();
      evalExpr->evaluate(evalMgr, coeffs);

    }

  
}

