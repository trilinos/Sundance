/* @HEADER@ */
/* @HEADER@ */


#include "SundanceTempStack.hpp"
#include "SundanceEvalVector.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;





TempStack::TempStack(int vecSize)
  : 
  vecSize_(vecSize),
  stack_(),
  numVecsAllocated_(0),
  numVecsAccessed_(0)
{}

TempStack::TempStack()
  : 
  vecSize_(0),
  stack_(),
  numVecsAllocated_(0),
  numVecsAccessed_(0)
{}


void TempStack::pushVectorData(const RefCountPtr<Array<double> >& vecData)
{
  stack_.push(vecData);
}

RefCountPtr<Array<double> > TempStack::popVectorData()
{
  RefCountPtr<Array<double> > data;
  if (stack_.empty())
    {
      numVecsAllocated_++;
      data = rcp(new Array<double>(vecSize_));
    }
  else
    {
      data = stack_.top();
      data->resize(vecSize_);
      stack_.pop();
    }
  numVecsAccessed_++;
  return data;
}



void TempStack::resetCounter()
{
  numVecsAllocated_=0;
  numVecsAccessed_=0;
}


