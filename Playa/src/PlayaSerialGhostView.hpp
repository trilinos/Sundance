/* @HEADER@ */
//
/* @HEADER@ */

#ifndef PLAYA_SERIAL_GHOSTVIEW_HPP
#define PLAYA_SERIAL_GHOSTVIEW_HPP

#include "PlayaDefs.hpp"
#include "PlayaGhostView.hpp"
#include "PlayaSerialVector.hpp"
#include "Teuchos_Utils.hpp"



namespace Playa
{
using namespace Teuchos;


/**
 * Dummy ghost element viewer for serial vectors. 
 */
class SerialGhostView : public GhostView<double>
{
public:
  /** */
  SerialGhostView(const RCP<SerialVector>& vec)
    : vec_(vec) 
    {;}

  /** virtual dtor */
  virtual ~SerialGhostView(){;}

  /** Indicate whether the given global index is accessible in this view */
  bool isAccessible(OrdType globalIndex) const 
    {return true;}

  /** get the element at the given global index */
  const double& getElement(OrdType globalIndex) const 
    {return (*vec_)[globalIndex];}

  /** get the batch of elements at the given global indices */
  void getElements(const OrdType* globalIndices, OrdType numElems,
    Array<double>& elems) const 
    {
      vec_->getElements(globalIndices, numElems, elems);
    }

  /** */
  void print(std::ostream& os) const {os << vec_->description();}
private:
  RCP<const SerialVector> vec_;
};
  
}

#endif
