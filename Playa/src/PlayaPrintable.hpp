/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_PRINTABLE_HPP
#define PLAYA_PRINTABLE_HPP

#include "PlayaDefs.hpp"
#include <iostream>

namespace Playa
{
  /**
   * Printable defines an interface for writing an object to a stream. 
   *
   * @author Kevin Long (kevin.long@ttu.edu)
   */
  class Printable
    {
    public:
      /** virtual dtor */
      virtual ~Printable() {;}

      /** abstract print function */
      virtual void print(std::ostream& os) const = 0 ;
    };
}


#endif
