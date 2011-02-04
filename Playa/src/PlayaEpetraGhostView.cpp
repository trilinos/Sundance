/* @HEADER@ */
//
/* @HEADER@ */

#include "PlayaEpetraGhostView.hpp"
#include "Epetra_Import.h"


namespace Playa
{

using namespace Teuchos;

const double& EpetraGhostView::getElement(OrdType globalIndex) const 
{
  const Epetra_BlockMap& myMap = ghostView_->Map();
  return (*ghostView_)[myMap.LID(globalIndex)];
}

void EpetraGhostView::getElements(const OrdType* globalIndices, int numElems,
  Array<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = ghostView_->Map();

  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*ghostView_)[myMap.LID(globalIndices[i])];
  }
}

void  EpetraGhostView::import(const Epetra_Import& importer,
  const Epetra_Vector& srcObject)
{
  /* If my vector does not yet exist, create it using the target map of the
   * importer */
  if (ghostView_.get()==0)
  {
    ghostView_ = rcp(new Epetra_Vector(importer.TargetMap()));
  }

  /* do the import */
  int ierr = ghostView_->Import(srcObject, importer, Insert);
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "ierr=" << ierr << " in EpetraGhostView::import()");
}

void EpetraGhostView::print(std::ostream& os) const
{
  if (ghostView_.get()==0) 
  {
    os << "[null Epetra ghost view]" << std::endl;
  }
  else
  {
    ghostView_->Print(os);
  }
}

}
