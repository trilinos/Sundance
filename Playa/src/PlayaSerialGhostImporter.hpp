/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_SERIAL_GHOSTIMPORTER_HPP
#define PLAYA_SERIAL_GHOSTIMPORTER_HPP

#include "PlayaDefs.hpp"
#include "PlayaGhostImporter.hpp"
#include "PlayaSerialGhostView.hpp"
#include "Teuchos_Utils.hpp"



namespace Playa
{
  using namespace Teuchos;


  /**
   * Ghost element importer for serial vectors. This class doesn't have
   * much to do, but is necessary to maintain a consistent interface.
   */
  class SerialGhostImporter : public GhostImporter<double>
    {
    public:
      /** */
      SerialGhostImporter(){;}
      /** virtual dtor */
      virtual ~SerialGhostImporter() {;}

      /** 
       * Import the ghost elements of the given vector
       * as specified during construction of this object. 
       */
      virtual void importView(const Vector<double>& x,
                              RCP<GhostView<double> >& ghostView) const ;

    private:
      
    };
  
}

#endif
