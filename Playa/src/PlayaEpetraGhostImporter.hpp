/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_EPETRAGHOSTIMPORTER_HPP
#define PLAYA_EPETRAGHOSTIMPORTER_HPP

#include "PlayaDefs.hpp"
#include "PlayaGhostImporter.hpp"
#include "PlayaGhostView.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Teuchos_Utils.hpp"



namespace Playa
{
  using namespace Teuchos;


  /**
   * Ghost element importer for Epetra vectors
   */
  class EpetraGhostImporter : public GhostImporter<double>
    {
    public:
      /** */
      EpetraGhostImporter(const RCP<const Epetra_Map>& epetraMap,
                          int nGhost,
                          const int* ghostElements);
      /** virtual dtor */
      virtual ~EpetraGhostImporter() {;}

      /** 
       * Import the ghost elements of the given vector
       * as specified during construction of this object. 
       */
      virtual void importView(const Vector<double>& x,
                              RCP<GhostView<double> >& ghostView) const ;

    private:
      RCP<const Epetra_Map> localMap_;

      RCP<const Epetra_Map> ghostMap_;

      RCP<Epetra_Import> importer_;
    };
  
}

#endif
