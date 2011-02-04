/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_GENERALIZED_INDEX_HPP
#define PLAYA_GENERALIZED_INDEX_HPP

#include "PlayaDefs.hpp"
#include <deque>

namespace Playa
{

/**
 * GeneralizedIndex is a locator for an element in 
 * vector with arbitrary block structure. Together with the processor's rank
 * it can uniquely specify an element in a distributed vector. 
 * This object is used in operations such as minloc and maxloc on
 * arbitrarily-structured vectors. 
 *
 * The implementation is a simple pair of (1) a deque giving the sequence of
 * block indices and (2) an integer giving the local offset within the
 * lowest-level block. The integer index is local, i.e., relative to the
 * block pointed to by the block index.
 */
class GeneralizedIndex
{
public:
  /** */
  GeneralizedIndex() : blockIndex_(), localIndex_(-1){}

  /** */
  void pushBlockIndex(int b) {blockIndex_.push_back(b);}

  /** */
  void popBlockIndex() {blockIndex_.pop_back();}

  /** */
  int readBlockIndex() const {return blockIndex_.back();}

  /** */
  void setLocalIndex(int i) {localIndex_=i;}

  /** */
  int readLocalIndex() const {return localIndex_;}

  /** */
  int depth() const {return blockIndex_.size();}

  /** */
  const std::deque<int>& blockIndex() const {return blockIndex_;}
  
  
private:
  std::deque<int> blockIndex_;
  int localIndex_;
};

}

#endif
