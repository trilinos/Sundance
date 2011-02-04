/* @HEADER@ */
//
/* @HEADER@ */

#include "PlayaEpetraGhostImporter.hpp"
#include "PlayaEpetraGhostView.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;


EpetraGhostImporter
::EpetraGhostImporter(const RCP<const Epetra_Map>& localMap,
  int nGhost,
  const int* ghostElements)
  : localMap_(localMap),
    ghostMap_(),
    importer_()
{
  if (false && nGhost==0)
  {
    ghostMap_ = localMap_;

    importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
  }
  else
  {
    //bvbw not used      int nGlobal = localMap_->NumGlobalElements();
    int nLocal = localMap_->NumMyElements();
    int nGhostView = nLocal+nGhost;
    std::vector<int> globalIndices(nGhostView);
    for (int i=0; i<nLocal; i++) globalIndices[i] = localMap_->GID(i);
    for (int i=0; i<nGhost; i++) globalIndices[i+nLocal] = ghostElements[i];

    const Epetra_Comm& comm = localMap_->Comm();

    ghostMap_ = rcp(new Epetra_Map(-1, nGhostView, 
        &(globalIndices[0]), 0, comm));

    importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
  }
}

void EpetraGhostImporter
::importView(const Vector<double>& x,
  RCP<GhostView<double> >& ghostView) const
{
  Tabs tab;

  /* If given an uninitialized ghost view, create a EpetraGhostView */
  if (ghostView.get()==0) 
  {
    ghostView = rcp(new EpetraGhostView());
  }

  /* Ensure that the ghost view contains an EpetraGhostView */
  EpetraGhostView* epgv 
    = dynamic_cast<EpetraGhostView*>(ghostView.get());

  TEST_FOR_EXCEPTION(epgv==0, std::runtime_error,
    "argument ghostView to EpetraGhostImporter::importView() "
    "could not be cast to a EpetraGhostView pointer");

  const Epetra_Vector& xVec = EpetraVector::getConcrete(x);

  Out::os() << tab << "EpGI importing" << std::endl;
  /* Do the import */
  epgv->import(*importer_, xVec);
}
    
}
