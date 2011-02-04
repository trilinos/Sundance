/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_OBJECT_WITH_VERBOSITY_HPP
#define PLAYA_OBJECT_WITH_VERBOSITY_HPP

#include "PlayaDefs.hpp"
#include <iostream>

namespace Playa
{
  /**
   * ObjectWithVerbosity provides a common interface for reading and setting
   * verbosity levels 
   *
   * @author Kevin Long (kevin.long@ttu.edu)
   */
  class ObjectWithVerbosity
    {
    public:
      /** */
      ObjectWithVerbosity(int verb=0) : verb_(verb) {}

      /** virtual dtor */
      virtual ~ObjectWithVerbosity() {;}

      /** Return the verbosity */
      virtual int verb() const {return verb_;}

      /** Set the verbosity */
      virtual void setVerb(int v) {verb_=v;}
    private:
      int verb_;
    };
}


#endif
