/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LOADABLEVECTOR_H
#define SUNDANCE_LOADABLEVECTOR_H


#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;

  namespace FrameworkInterface
    {
      /**
       * LoadableVector is the interface through which a framework
       * inserts values into a symbolic expression.
       */
      class LoadableVector
        {
        public:
          /** */
          LoadableVector();

          /** */
          virtual ~LoadableVector(){;}

          /** Change the size of the vector to newSize */
          virtual void resize(int newSize) = 0 ;

          /** Return the length of the vector */
          virtual int length() const = 0 ;

          /** Set the i-th element to x */
          virtual void setElement(int i, const double& x) = 0 ;

          /** Return a pointer to the physical start of the vector. */
          virtual double* const start() = 0 ;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
