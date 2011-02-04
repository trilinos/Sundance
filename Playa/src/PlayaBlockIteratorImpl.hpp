/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_BLOCKITERATORIMPL_HPP
#define PLAYA_BLOCKITERATORIMPL_HPP


#include "PlayaBlockIteratorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaExceptions.hpp"

namespace Playa
{

template <class Scalar> inline
bool BlockIterator<Scalar>
::operator==(const BlockIterator<Scalar>& other) const
{
  if (debug())
  {
    Out::os() << "comparing: LHS=" << *this
              << ", RHS=" << other << std::endl;
  }
  if (this->atEnd_ && other.atEnd_) return true;
  if (this->atEnd_ != other.atEnd_) return false;
  if (this->index_.size() != other.index_.size()) return false;

  for (unsigned int i=0; i<this->index_.size(); i++)
  {
    if (this->index_[i] != other.index_[i]) return false;
  }

  return true;
}

template <class Scalar> inline
bool BlockIterator<Scalar>
::operator<(const BlockIterator<Scalar>& other) const
{
  if (debug())
  {
    Out::os() << "comparing (<): LHS=" << *this
              << ", RHS=" << other << std::endl;
  }
  TEST_FOR_EXCEPTION(this->space() != other.space(),
    RuntimeError, "Attempt to compare block iterators attached to "
    "two different spaces");

  if (!this->atEnd_ && other.atEnd_) return true;
  if (this->atEnd_ && !other.atEnd_) return false;
  if (this->atEnd_ && other.atEnd_) return false;

  int d = std::min(this->index_.size(), other.index_.size());
  

  for (int i=0; i<d; i++)
  {
    if (this->index_[i] < other.index_[i]) return true;
    if (this->index_[i] > other.index_[i]) return false;
  }

  return false;
}



template <class Scalar> inline
BlockIterator<Scalar> BlockIterator<Scalar>::operator++(int)
{
  if (debug())
  {
    Out::os() << "iter++" << std::endl;
  }
  BlockIterator<Scalar> old = *this;

  TEST_FOR_EXCEPTION(this->atEnd_, RuntimeError,
    "attempt to advance a BlockIterator beyond end");

  int depth = this->index_.size();
  TEST_FOR_EXCEPTION(depth <= 0, RuntimeError, 
    "empty index stack in BlockIterator");

  atEnd_ = !advance(depth-1);

  if (debug())
  {
    Out::os() << "old=" << old << ", new=" << *this << std::endl;
  }
  return old;
}

template <class Scalar> inline
bool BlockIterator<Scalar>::advance(int level) 
{
  if (level < 0) return false; 

  this->index_[level]++;
  std::deque<int> base = this->index_;
  base.pop_back();
  /* If we're at the end of this level, drop back a level */
  if (this->index_[level] >= this->space().getBlock(base).numBlocks())
  {
    this->index_.pop_back();
    return this->advance(level-1);
  }
  else /* go to the start of the next block */
  {
    goToStart(this->space().getBlock(*this), this->index_);
  }
  return true;
}

template <class Scalar> inline
void BlockIterator<Scalar>
::goToStart(const VectorSpace<Scalar>& space,
  std::deque<int>& pos) const
{
  pos.push_back(0);
  if (space.isBlockSpace())
  {
    goToStart(space.getBlock(0), pos);
  }
}


template <class Scalar> inline
const VectorSpace<Scalar>& BlockIterator<Scalar>::space() const
{
  return *(this->space_);
}


template <class Scalar> inline
BlockIterator<Scalar>::BlockIterator(
  const VectorSpace<Scalar>& space,
  bool atEnd
  )
  : space_(rcp(new VectorSpace<Scalar>(space))),
    index_(),
    atEnd_(atEnd)
{
  if (!atEnd) 
  {
    goToStart(space, index_);
  }
}


template <class Scalar> inline
std::ostream& BlockIterator<Scalar>::print(std::ostream& os) const 
{
  os << "BlockIterator";
  if (this->atEnd_) 
  {
    os << "[end]";
  }
  else
  {
    os << "[index={";
    for (unsigned int i=0; i<index_.size(); i++) 
    {
      if (i>0U) os << ", ";
      os << index_[i];
    }
    os << "}]";
  }
  return os;
}



}





#endif
