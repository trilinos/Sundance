/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_BLOCKVECTORSPACEBASEDECL_HPP
#define PLAYA_BLOCKVECTORSPACEBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"

namespace Playa
{

/**
 * This is a base class for a blocked vector space. It assumes nothing about
 * the physical storage of the blocks. 
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class BlockVectorSpaceBase : public VectorSpaceBase<Scalar>
{
public:
  /** virtual dtor */
  virtual ~BlockVectorSpaceBase() {;}

  /** Compute dimension. The default implementation sums the dimensions
   * of all blocks. */
  virtual int dim() const ;

  /** Count the locally owned elements */
  virtual int numLocalElements() const ;

  /** Check compatibility with another space. The default implementation
   * loops over blocks checking compatiblilty at each index. */
  virtual bool isCompatible(const VectorSpaceBase<Scalar>* other) const ;

  /** Get a block specified by an integer index. This function
   * should hrow an exception if the index is out of range 
   */
  virtual const VectorSpace<Scalar>& getBlock(int b) const = 0 ; 

  /** Write a description by recursivle describing the blocks */
  virtual std::string description() const ;

};






}

#endif
