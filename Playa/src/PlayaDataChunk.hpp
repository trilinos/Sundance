/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_DATACHUNK_HPP
#define PLAYA_DATACHUNK_HPP

#include "PlayaDefs.hpp"

namespace Playa
{
/**
 * This class wraps a pointer to a single contiguous chunk of scalar data 
 * 
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class ConstDataChunk
{
public:
  /** */
  ConstDataChunk(int size, const Scalar* values)
    : size_(size), values_(values) {}

  /** */
  const Scalar* values() const {return values_;}

  /** */
  int size() const {return size_;}

private:
  int size_;
  const Scalar* values_;
};

/**
 * This class wraps a pointer to a single contiguous chunk of scalar data 
 * 
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class NonConstDataChunk
{
public:
  /** */
  NonConstDataChunk(int size, Scalar* values)
    : size_(size), values_(values) {}

  /** */
  const Scalar* values() const {return values_;}

  /** */
  Scalar* values() {return values_;}

  /** */
  int size() const {return size_;}
private:
  int size_;
  Scalar* values_;
};

  
  
}

#endif
