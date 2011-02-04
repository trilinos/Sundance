/* @HEADER@ */

 /* @HEADER@ */

#ifndef PLAYA_BLOCKOPERATORBASE_DECL_HPP
#define PLAYA_BLOCKOPERATORBASE_DECL_HPP

#include "PlayaDefs.hpp"


namespace Playa
{


template <class Scalar>
class LinearOperator;

/**
 * Class BlockOperatorBase provides an abstract interface for accessing
 * blocks of block operators
 *
 * @author Paul T Boggs (ptboggs@sandia.gov)
 */
template <class Scalar>
class BlockOperatorBase 
{
public:

  /** */
  virtual int numBlockRows() const = 0 ;

  /** */
  virtual int numBlockCols() const = 0 ;

  /** */
  virtual const LinearOperator<Scalar>& getBlock(int i, int j) const = 0 ;
  /** */
  virtual LinearOperator<Scalar> getNonconstBlock(int i, int j) = 0 ;

}; 

/**
 * Class SetableBlockOperatorBase provides an abstract interface for setting
 * blocks of block operators
 *
 * @author Paul T Boggs (ptboggs@sandia.gov)
 */
template <class Scalar>
class SetableBlockOperatorBase : public BlockOperatorBase<Scalar>
{
public:

  /** */
  virtual void setBlock(int i, int j, const LinearOperator<Scalar>& Aij) = 0 ;

  /** Callback to finalize block fill. The default is a no-op. */
  virtual void endBlockFill() {;}

  

}; 





}

#endif
