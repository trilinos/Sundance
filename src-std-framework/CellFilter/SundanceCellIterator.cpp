/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellIterator.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

CellIterator::CellIterator()
  :  isImplicit_(true),
    currentLID_(-1),
    reorderer_(0),
    iter_()
{;}

CellIterator::CellIterator(const Mesh& mesh, 
                           int cellDim, 
                           CellIteratorPos pos)
  : isImplicit_(true),
    currentLID_(-1),
    reorderer_(mesh.reorderer()),
    iter_()
{
  if (cellDim == mesh.spatialDim() && reorderer_ != 0)
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = reorderer_->begin(); 
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(verbosity<CellFilter>() > VerbMedium, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }
  else 
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = 0;
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(verbosity<CellFilter>() > VerbMedium, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }


}



CellIterator::CellIterator(const Set<int>* cells, CellIteratorPos pos)
  : isImplicit_(false),
    currentLID_(-1),
    reorderer_(0),
    iter_()
{
  switch(pos)
    {
    case Begin:
      iter_ = cells->begin();
      break;
    case End:
      iter_ = cells->end();
    }
  SUNDANCE_OUT(verbosity<CellFilter>() > VerbMedium, 
               "created explicit cell iterator");
}



    
