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


static Time& getLocalValsTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("DF getLocalValues()"); 
  return *rtn;
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const double& constantValue,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{
  vector_.setToConstant(constantValue);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Vector<double>& vector,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(vector),
    ghostView_(),
    ghostsAreValid_(false)
{}

void DiscreteFunction::setVector(const Vector<double>& vec) 
{
  ghostsAreValid_ = false;
  vector_ = vec;
}

void DiscreteFunction::updateGhosts() const
{
  if (!ghostsAreValid_)
    {
      space_.importGhosts(vector_, ghostView_);
      ghostsAreValid_ = true;
    }
}


void DiscreteFunction::getLocalValues(int cellDim, 
                        const Array<int>& cellLID,
                        Array<double>& localValues) const 
{
  TimeMonitor timer(getLocalValsTimer());

  updateGhosts();

  const RefCountPtr<DOFMapBase>& map = space_.map();
  static Array<int> dofs;
  static Array<int> indices;
  int nNodes;
  map->getDOFsForCellBatch(cellDim, cellLID, dofs, nNodes);
  int nFunc = space_.nFunc();
  int nCells = cellLID.size();
  localValues.resize(nFunc*cellLID.size()*nNodes);
  indices.resize(dofs.size());
  

  for (int c=0; c<cellLID.size(); c++)
    {
      for (int f=0; f<nFunc; f++)
        {
          for (int n=0; n<nNodes; n++)
            {
              indices[c*nFunc*nNodes + nFunc*n + f]
                = dofs[(f*nCells + c)*nNodes + n];
            }
        }
    }
  ghostView_->getElements(&(indices[0]), indices.size(), localValues);
}


const DiscreteFunction* DiscreteFunction::discFunc(const Expr& expr)
{
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



DiscreteFunction* DiscreteFunction::discFunc(Expr& expr)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



