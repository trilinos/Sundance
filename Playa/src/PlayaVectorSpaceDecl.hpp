/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORSPACEDECL_HPP
#define PLAYA_VECTORSPACEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"
#include "PlayaBlockIteratorDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 *  User-level VectorSpace class.
 */
template <class Scalar>
class VectorSpace : public Playa::Handle< const VectorSpaceBase<Scalar> >
{
public:
  HANDLE_CTORS(VectorSpace<Scalar>, const VectorSpaceBase<Scalar>);
    
  /** Create a new element of this vector space */
  Vector<Scalar>  createMember() const ;

  /** Return the dimension of the space */
  int dim() const {return this->ptr()->dim();}

  /** Return the lowest global index accessible on this processor */
  int baseGlobalNaturalIndex() const ;

  /** Return the number of elements owned by this processor */
  int numLocalElements() const ;

  /** Return the MPI communicator */
  const MPIComm& comm() const {return this->ptr()->comm();}

  /** Check compatibility with another space. */
  bool isCompatible(const VectorSpace<Scalar>& vecSpc) const; 


  /** test equality between two spaces */
  bool operator==(const VectorSpace<Scalar>& other) const ;


  /** test inequality of two spaces */
  bool operator!=(const VectorSpace<Scalar>& other) const ;


  /** test whether the space contains a given vector */
  bool contains(const Vector<Scalar>& vec) const ;


  /** return the number of subblocks at the highest level. */
  int numBlocks() const ;

  /** indicate whether I am a block vector space */
  bool isBlockSpace() const ;

  /** get the i-th subblock */
  const VectorSpace<Scalar>& getBlock(int i) const ;

  /** get a subblock as specified by a block iterator */
  const VectorSpace<Scalar>& getBlock(const BlockIterator<Scalar>& iter) const ;

  /** get a subblock as specified by a deque of indices */
  const VectorSpace<Scalar>& getBlock(const std::deque<int>& iter) const ;

  /** */
  BlockIterator<Scalar> beginBlock() const ;

  /** */
  BlockIterator<Scalar> endBlock() const ;

  /** */
  int mapToGNI(const BlockIterator<Scalar>& b, int indexWithinBlock) const ;

  /** */
  bool containsGNI(int gni) const ;

  /** */
  void getBlockAndOffsetFromGNI(int gni,
    BlockIterator<Scalar>& block, int& indexWithinBlock) const ;
  
protected:
  

};

template <class Scalar>
STREAM_OUT(VectorSpace<Scalar>)

}


#endif
