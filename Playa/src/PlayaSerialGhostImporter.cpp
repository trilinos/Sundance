/* @HEADER@ */
//
 /* @HEADER@ */

#include "PlayaSerialGhostImporter.hpp"
#include "PlayaSerialGhostView.hpp"
#include "PlayaSerialVector.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif

using namespace Teuchos;
using namespace Playa;



void SerialGhostImporter
::importView(const Vector<double>& x,
             RCP<GhostView<double> >& ghostView) const
{
  RCP<SerialVector> xPtr = rcp_dynamic_cast<SerialVector>(x.ptr());
  ghostView = rcp(new SerialGhostView(xPtr));
}
    
