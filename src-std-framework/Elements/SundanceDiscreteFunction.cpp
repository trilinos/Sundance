/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFunction.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"

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
                                   const double& constantValue,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector())
{
  vector_.setToConstant(constantValue);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Vector<double>& vector,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(vector)
{}


void DiscreteFunction::getLocalValues(int cellDim, 
                        const Array<int>& cellLID,
                        Array<double>& localValues) const 
{
  const RefCountPtr<DOFMapBase>& map = space_.map();
  Array<int> dofs;
  map->getDOFsForCell(cellDim, 0, 0, dofs);
  int nNodes = dofs.size();
  int nFunc = space_.nFunc();
  localValues.resize(nFunc*cellLID.size()*nNodes);

  for (int c=0; c<cellLID.size(); c++)
    {
      for (int f=0; f<nFunc; f++)
        {
          map->getDOFsForCell(cellDim, cellLID[c], f, dofs);

          for (int n=0; n<nNodes; n++)
            {
              localValues[c*nFunc*nNodes + nFunc*n + f] 
                = vector_.getElement(dofs[n]);
            }
        }
    }
}

// DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
//                                    const Expr& expr,
//                                    const string& name)
//   : DiscreteFunctionStub(name, space.nFunc()), 
//     FuncWithBasis(space.basis()),
//     space_(space),
//     vector_(space_.createVector())
// {
//  //  Tabs tab;

//   RefCountPtr<Array<int> > workSet = rcp(new Array<int>());
//   workSet->reserve(workSetSize());
//   RefCountPtr<Array<double> > localValues = rcp(new Array<double>());

//   TSFExtended::LoadableVector<double>* vec 
//     = dynamic_cast<TSFExtended::LoadableVector<double>* >(vector_.ptr().get());

//   TEST_FOR_EXCEPTION(vec==0, RuntimeError,
//                      "vector is not loadable in DiscreteFunction ctor");


//   CellFilter filter = new MaximalCellFilter();

//   CellSet cells = filter.getCells(space_.mesh());
//   int cellDim = filter.dimension(mesh_);
//   CellType cellType = space_.mesh().cellType(cellDim);

//   CellIterator iter=cells.begin();
//   int workSetCounter = 0;

//   while (iter != cells.end())
//     {
//       Tabs tab1;
//       /* build up the work set */
//       workSet->resize(0);
//       for (int c=0; c<workSetSize() && iter != cells.end(); c++, iter++)
//         {
//           workSet->append(*iter);
//         }
//       SUNDANCE_OUT(verbosity() > VerbMedium,
//                    tab1 << "doing work set " << workSetCounter
//                    << " consisting of " << workSet->size() << " cells");
//       workSetCounter++;

//       mediator->setCellBatch(workSet);

//       evalExpr->flushResultCache();
//       evalExpr->evaluate(evalMgr, coeffs);

//     }

  
// }

