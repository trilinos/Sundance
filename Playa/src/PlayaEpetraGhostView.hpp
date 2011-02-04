/* @HEADER@ */
//
/* @HEADER@ */

#ifndef PLAYA_EPETRA_GHOST_VIEW_HPP
#define PLAYA_EPETRA_GHOST_VIEW_HPP

#include "PlayaDefs.hpp"
#include "PlayaGhostImporter.hpp"
#include "PlayaGhostView.hpp"
#include "Epetra_Vector.h"
#include "Teuchos_Utils.hpp"



namespace Playa
{
using namespace Teuchos;


/**
 * Ghost element viewer for Epetra vectors
 */
class EpetraGhostView : public GhostView<double>
{
public:
  /** */
  EpetraGhostView()
    : ghostView_() 
    {;}

  /** virtual dtor */
  virtual ~EpetraGhostView(){;}

  /** Indicate whether the given global index is accessible in this view */
  bool isAccessible(OrdType globalIndex) const 
    {return ghostView_->Map().MyGID(globalIndex);}

  /** get the element at the given global index */
  const double& getElement(OrdType globalIndex) const ;

  /** get the batch of elements at the given global indices */
  void getElements(const OrdType* globalIndices, int numElems,
    Array<double>& elems) const ;

  /** */
  void import(const Epetra_Import& importer,
    const Epetra_Vector& srcObject);

  /** */
  void print(std::ostream& os) const ;
private:
  RCP<Epetra_Vector> ghostView_;
};
  
}

#endif
