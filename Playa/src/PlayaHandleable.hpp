/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_HANDLEABLE_HPP
#define PLAYA_HANDLEABLE_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"

#define GET_RCP(Base) \
  /** Handleable<##Base> interface */ \
  virtual Teuchos::RCP<Base > getRcp() {return rcp(this);}

namespace Playa
{
using namespace Teuchos;

/**
 * Class Handleable provides an abstract interface for polymorphic
 * conversion from raw pointers to smart pointers. Recall from the
 * Teuchos RefCountPtr documentation that one should never create
 * directly a smart pointer from a raw pointer; rather, smart pointers
 * should be created through a call to rcp(). The type of the argument
 * to rcp() must be known at compile time. This makes the syntax
 * \code
 * Handle h = new Derived();
 * \endcode
 * impossible with the straightforward implementation in which Handle takes
 * a raw pointer to a Base. In order to preserve this clean syntax, we
 * require any handles supporting this syntax to take a raw
 * pointer to a Handleable<Base>, where Handleable<Base> provides a 
 * getRcp() method which returns the result of a call to rcp() on this.
 */
template <class Base>
class Handleable
{
public:
  /** Virtual dtor */
  virtual ~Handleable(){;}

  /** Return a safely-created RefCountPtr to the base type */
  virtual RCP<Base> getRcp() = 0 ;
    
};
  
}




#endif
