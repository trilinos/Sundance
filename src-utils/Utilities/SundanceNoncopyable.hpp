/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NONCOPYABLE_H
#define SUNDANCE_NONCOPYABLE_H

#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  /**
   * Noncopyable has a private copy ctor and assignment operator, which
   * prevents any classes dervied from it from being copied. 
   *
   * Based on the Noncopyable class from www.boost.org.
   */
  class Noncopyable
  {
  public:
    Noncopyable(){;}

  private:
    /** private copy ctor */
    Noncopyable(const Noncopyable& other);

    /** private assignment operator */
    Noncopyable& operator=(const Noncopyable& other);
  };
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
