/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_MATRIXFACTORY_HPP
#define PLAYA_MATRIXFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaLinearOperatorDecl.hpp"

namespace Playa
{
  /** 
   * MatrixFactory is an abstract builder for empty matrices
   */
  template <class Scalar>
  class MatrixFactory
  {
  public:
    /** Virtual dtor */
    virtual ~MatrixFactory(){;}

    /** */
    virtual LinearOperator<Scalar> createMatrix() const = 0 ;
  };
}

#endif
