/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_DENSE_SERIAL_MATRIXFACTORY_HPP
#define PLAYA_DENSE_SERIAL_MATRIXFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaMatrixFactory.hpp"
#include "PlayaSerialVectorSpace.hpp"

namespace Playa
{
  /** 
   * 
   */
class DenseSerialMatrixFactory : public MatrixFactory<double>
{
public:
  /** */
  DenseSerialMatrixFactory(
    const VectorSpace<double>& domain,
    const VectorSpace<double>& range);

  /** Virtual dtor */
  virtual ~DenseSerialMatrixFactory(){;}
  
  /** */
  virtual LinearOperator<double> createMatrix() const ;

public:
  VectorSpace<double> domain_;
  VectorSpace<double> range_;
};
}

#endif
